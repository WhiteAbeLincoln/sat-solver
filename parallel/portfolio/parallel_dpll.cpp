/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::*Full%20Source][Full Source:1]] */
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <stack>
#include <tuple>
#include <sstream>
#include <unistd.h>
#include <string>
#include <cmath>
#include <cstdlib>
#include <random>
#include <mpi.h>
#define MCW MPI_COMM_WORLD

int sign(int x) {
  return ( (x > 0) ? 1
         : (x < 0) ? -1
         : 0);
}
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::read_input][read_input]] */
auto read_input() {
  std::vector<int> f;
  for (std::string l; std::getline(std::cin, l);) {
    if (l.empty()) continue;
    std::stringstream ss(l);
    std::string word;
    ss >> word;
    if (word == "c") continue;
    if (word == "p") {
      ss >> word;
      if (word != "cnf") throw std::invalid_argument("Data must be in cnf format, got " + word);
      continue;
    }
    do {
      if (word == "%") return f;
      int v = std::stoi(word);
      f.push_back(v);
    } while (ss >> word);
  }
  
  return f;
}
/* read_input ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::branch_enum_def][branch_enum_def]] */
enum class BranchRule { dlis, dlcs, jw, jw2, dsj, rand };
/* branch_enum_def ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::formula][formula]] */
  struct Formula {
    struct ClauseData {
      int n_t = 0;
      int n_f = 0;
      std::vector<int> literals;
      int orig_len;
      bool sat() { return n_t >= 1; }
      bool unsat() { return n_f == orig_len; }
      bool unit() { return n_t == 0 && n_f == (orig_len - 1); }
    };
    struct LitData {
      int assn = -1;
      std::vector<int> pos_clauses;
      std::vector<int> neg_clauses;
      bool pure() {
        return assn == -1 && (pos_clauses.size() == 0 || neg_clauses.size() == 0);
      }
    };
    std::vector<ClauseData> clauses;
    std::unordered_map<int, LitData> literals;
    int remaining;
    void add_literal(int, int);
    Formula(std::vector<int>); 
  };
/* formula ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::remove_satisfied][remove_satisfied]] */
void remove_satisfied(Formula& f, int d) {
  auto& clause = f.clauses[d];
  clause.n_t++;
  f.remaining--;
  auto lits = clause.literals;
  for (auto l : lits) {
    auto s = sign(l);
    auto& lit = f.literals[l*s];
    if (s == 1) {
      auto& p = lit.pos_clauses;
      p.erase(std::remove(p.begin(), p.end(), d), p.end());
    } else {
      auto& n = lit.neg_clauses;
      n.erase(std::remove(n.begin(), n.end(), d), n.end());
    }
  }
  clause.literals.clear();
}
/* remove_satisfied ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::set_var][set_var]] */
bool set_var(Formula& f, int l) {
  auto s = sign(l);
  auto pos = l*s;
  auto& lit = f.literals[pos];
  if (lit.assn != -1) throw std::runtime_error("literal already assigned");
  lit.assn = (s == 1) ? 1 : 0;
  auto sat_c = (lit.assn == 1) ? lit.pos_clauses : lit.neg_clauses;
  auto& unsat_c = (lit.assn == 0) ? lit.pos_clauses : lit.neg_clauses;
  for (auto cidx : sat_c) remove_satisfied(f, cidx);
  for (auto cidx : unsat_c) {
    auto& clause = f.clauses[cidx];
    clause.n_f++;
    clause.literals.erase(std::remove(clause.literals.begin(),
                                      clause.literals.end(),
                                      (lit.assn == 0) ? pos : -pos),
                          clause.literals.end());
    if (clause.literals.size() == 0) return false;
  }
  unsat_c.clear();
  return true;
}
/* set_var ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::pure_literal][pure_literal]] */
void pure_literal_assign(Formula& f, Formula::LitData& data) {
  auto pos_size = data.pos_clauses.size();
  auto s = (pos_size == 0) ? -1 : 1;
  auto lclauses = (s == 1) ? data.pos_clauses : data.neg_clauses;
  data.assn = (s == 1) ? 1 : 0;
  for (auto cidx : lclauses) remove_satisfied(f, cidx);
}
/* pure_literal ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::unit_propogate][unit_propogate]] */
bool unit_propogate(Formula& f, Formula::ClauseData clause) {
  return set_var(f, clause.literals[0]);
}
/* unit_propogate ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::branch_helper][branch_helper]] */
int nclauses(Formula f, int k, int u) {
  auto s = sign(u);
  auto& lit = f.literals[u*s];
  auto cs = (s == 1) ? lit.pos_clauses : lit.neg_clauses;
  int counter = 0;
  if (k == -1) return cs.size();
  for (auto c : cs) {
    if (f.clauses[c].literals.size() == (unsigned int)k) counter++;
  }
  return counter;
}
int get_largest_k(Formula f) {
  return std::max_element(f.clauses.begin(), f.clauses.end(),
                  [](auto a, auto b) {
                    return a.literals.size() < b.literals.size();
                  })->literals.size();
}
int apply_rule(Formula f, std::function<std::tuple<int,int,int>(Formula, int)> rule) {
  int maximum = 0;
  int curr = 0;
  for (auto l : f.literals) {
    if (l.second.assn != -1) continue;
    auto [wp, wn, phi] = rule(f, l.first);
    if (phi >= maximum) {
      curr = wp >= wn ? l.first : -l.first;
      maximum = phi;
    }
  }
  return curr;
}
/* branch_helper ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::dlis][dlis]] */
auto dlis(Formula f, int l) {
  int wp = nclauses(f, -1, l);
  int wn = nclauses(f, -1, -l);
  return std::make_tuple(wp, wn, std::max(wp, wn));
}
/* dlis ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::dlcs][dlcs]] */
auto dlcs(Formula f, int l) {
  int wp = nclauses(f, -1, l);
  int wn = nclauses(f, -1, -l);
  return std::make_tuple(wp, wn, wp + wn);
}
/* dlcs ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::jw][jw]] */
auto jw(Formula f, int l) {
  auto largest_k = get_largest_k(f);
  int wp = 0;
  int wn = 0;
  for (int k = 1; k <= largest_k; ++k) {
    wp += std::pow(2, -k) * nclauses(f, k, l);
    wn += std::pow(2, -k) * nclauses(f, k, -l);
  }
  return std::make_tuple(wp, wn, std::max(wp, wn));
}
/* jw ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::jw2][jw2]] */
auto jw2(Formula f, int l) {
  auto largest_k = get_largest_k(f);
  int wp = 0;
  int wn = 0;
  for (int k = 1; k <= largest_k; ++k) {
    wp += std::pow(2, -k) * nclauses(f, k, l);
    wn += std::pow(2, -k) * nclauses(f, k, -l);
  }
  return std::make_tuple(wp, wn, wp + wn);
}
/* jw2 ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::dsj][dsj]] */
auto dsj(Formula f, int l) {
  auto largest_k = get_largest_k(f);
  int wp = 4*nclauses(f, 2, l) + 2*nclauses(f, 3, l);
  int wn = 4*nclauses(f, 2, -l) + 2*nclauses(f, 3, -l);
  for (int k = 4; k <= largest_k; ++k) {
    wp += nclauses(f, k, l);
    wn += nclauses(f, k, -l);
  }
  return std::make_tuple(wp, wn, (wp+1)*(wn+1));
}
/* dsj ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::get_branching][get_branching]] */
std::string branch_rule_name(BranchRule rule) {
  switch (rule) {
    case BranchRule::dlis:
      return "dlis";
    case BranchRule::dlcs:
      return "dlcs";
    case BranchRule::jw:
      return "jw";
    case BranchRule::jw2:
      return "jw2";
    case BranchRule::dsj:
      return "dsj";
    case BranchRule::rand:
      return "rand";
  }
  throw std::runtime_error("branch_rule_name didn't handle all cases");
}
int get_branching_variable(Formula f, BranchRule rule) {
  int curr = 0;
  switch (rule) {
    case BranchRule::dlis:
      curr = apply_rule(f, &dlis);
      break;
    case BranchRule::dlcs:
      curr = apply_rule(f, &dlcs);
      break;
    case BranchRule::jw:
      curr = apply_rule(f, &jw);
      break;
    case BranchRule::jw2:
      curr = apply_rule(f, &jw2);
      break;
    case BranchRule::dsj:
      curr = apply_rule(f, &dsj);
      break;
    case BranchRule::rand: {
      bool positive = std::rand() % 2;
      int i = 0;
      long unsigned int iter = 0;
      do {
        i = 1 + std::rand() % f.literals.size();
        iter++;
      } while (f.literals[i].assn != -1 && iter < f.literals.size());
      if (i == 0) throw std::runtime_error("random branch failed");
      curr = (positive ? i : -i);
      break;
    }
    default:
      throw std::runtime_error("get_branching_variable didn't handle all cases");
  }
  if (curr == 0)
    throw std::runtime_error("branching heuristic failed: " + branch_rule_name(rule));
  return curr;
}
BranchRule branch_rule_int(int i) {
  return (i > static_cast<int>(BranchRule::rand)) ? BranchRule::rand : static_cast<BranchRule>(i);
}
/* get_branching ends here */
/* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::dpll][dpll]] */
  std::tuple<bool, Formula> dpll(Formula& f, BranchRule rule, std::function<bool()> terminate) {
    if (terminate()) return {false, f};
    
    for (auto&& [_f, l] : f.literals)
      if (l.pure()) pure_literal_assign(f, l);

    for (auto& c : f.clauses) {
      if (c.sat()) continue;
      if (c.literals.size() == 0) return {false, f};
      if (c.unit())
        if (!unit_propogate(f, c)) return {false, f};
    }

    if (f.remaining == 0) return {true, f};

    auto l = get_branching_variable(f, rule);
    Formula oldf(f);
    set_var(f, l);
    auto [res, ff] = dpll(f, rule, terminate);
    if (res) return {res, ff};

    f = oldf;
    set_var(f, -l);
    return dpll(f, rule, terminate);
  }
/* dpll ends here */
void Formula::add_literal(int l, int cpos) {
  auto s = sign(l);
  auto pos = l*s;
  auto it = this->literals.find(pos);
  LitData data;
  if (it != this->literals.end()) data = it->second;
  if (s == 1) {
    data.pos_clauses.push_back(cpos);
  } else {
    data.neg_clauses.push_back(cpos);
  }
  this->literals.insert_or_assign(pos, data);
}
Formula::Formula(std::vector<int> formula) {
  ClauseData cd;
  for (auto l : formula) {
    if (l == 0) {
      cd.orig_len = cd.literals.size();
      this->clauses.push_back(cd);
      cd = ClauseData{};
    } else {
      cd.literals.push_back(l);
      auto cpos = this->clauses.size();
      this->add_literal(l, cpos);
    }
  }
  this->remaining = this->clauses.size();
}

int main(int argc, char** argv) {
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MCW, &rank);
  MPI_Comm_size(MCW, &size);

  /* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::data_sharing][data_sharing]] */
  int* form;
  int form_c;
  std::vector<int> f;
  if (rank == 0) {
    f = read_input();
    form = f.data();
    form_c = f.size();
  }
  MPI_Bcast(&form_c, 1, MPI_INT, 0, MCW);
  if (rank != 0) form = (int*)malloc(sizeof(int) * form_c);
  MPI_Bcast(form, form_c, MPI_INT, 0, MCW);
  if (rank != 0) f = std::vector<int>(form, form + form_c);
  Formula formula(f);
  /* data_sharing ends here */
  /* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::worker_process][worker_process]] */
  if (rank != 0) {
    std::srand(rank);
    auto rule = branch_rule_int(rank - 1);
    int early_term = 0;
    auto [sat, finalf] = dpll(formula, rule,
                              [&]() {
                                MPI_Iprobe(0, 0, MCW, &early_term, MPI_STATUS_IGNORE);
                                return static_cast<bool>(early_term);
                              });
  
    if (sat) {
      std::vector<int> assn;
      for (auto l : finalf.literals) {
        if (l.second.assn == 1) assn.push_back(l.first);
      }
      assn.push_back(0);
      for (auto l : finalf.literals) {
        if (l.second.assn == 0) assn.push_back(l.first);
      }
  
      MPI_Send(assn.data(), assn.size(), MPI_INT, 0, 0, MCW);
    } else if (!early_term) {
      int data = 0;
      MPI_Send(&data, 1, MPI_INT, 0, 0, MCW);
    }
  }
  /* worker_process ends here */
  /* [[file:~/Documents/Projects/sat-solver/parallel/portfolio/portfolio.org::master_process][master_process]] */
  if (rank == 0) {
    int* buff = (int*)malloc(sizeof(int) * (formula.literals.size() + 1));
    MPI_Status status;
    MPI_Recv(buff, formula.literals.size() + 1, MPI_INT, MPI_ANY_SOURCE, 0, MCW, &status);
    int solver_rank = status.MPI_SOURCE;
    for (int i = 1; i < size; ++i) {
      if (i == solver_rank) continue;
      int data = 0;
      MPI_Send(&data, 1, MPI_INT, i, 0, MCW);
    }
    
    if (buff[0] == 0) {
      std::cout << "Formula is: UNSAT" << std::endl;
      std::cout << "Solved By process " << solver_rank << " with branching strategy "
                << branch_rule_name(branch_rule_int(solver_rank-1)) << std::endl;
    } else {
      std::cout << "Formula is: SAT" << std::endl;
      std::cout << "Solved By process " << solver_rank << " with branching strategy "
                << branch_rule_name(branch_rule_int(solver_rank-1)) << std::endl;
      std::cout << "Variables assigned TRUE:" << std::endl;
      for (long unsigned int i = 0; i < formula.literals.size() + 1; ++i) {
        if (buff[i] == 0)
          std::cout << "\nVariables assigned FALSE:" << std::endl;
        else
          std::cout << buff[i] << " ";
      }
    }
  }
  /* master_process ends here */
  
  MPI_Finalize();
  return 0;
}
/* Full Source:1 ends here */
