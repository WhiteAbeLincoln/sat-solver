/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::*Full%20Source][Full Source:1]] */
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

bool VERBOSE = false;

/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::sign][sign]] */
int sign(int x) {
  return ( (x > 0) ? 1
         : (x < 0) ? -1
         : 0);
}
/* sign ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::read_input][read_input]] */
auto read_input() {
  bool cnf_mode = false;
  std::vector<std::vector<int>> f;
  std::vector<int> clause;
  for (std::string l; std::getline(std::cin, l);) {
    if (l.empty()) continue;
    std::stringstream ss(l);
    std::string word;
    ss >> word;
    if (word == "c") continue;
    if (word == "p") {
      ss >> word;
      if (word != "cnf") throw std::invalid_argument("Data must be in cnf format, got " + word);
      cnf_mode = true;
      continue;
    }
    do {
      if (word == "%") return f;
      int v = std::stoi(word);
      if (v == 0) {
        f.push_back(clause);
        clause.clear();
      } else {
        clause.push_back(v);
      }
    } while (ss >> word);
    if (!cnf_mode) { 
      f.push_back(clause);
      clause.clear();
    }
  }
  
  return f;
}
/* read_input ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::branch_enum_def][branch_enum_def]] */
enum class BranchRule { dlis, dlcs, jw, jw2, dsj };
/* branch_enum_def ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::formula][formula]] */
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
  void add_literal(int l, int cpos) {
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
  Formula(std::vector<std::vector<int>> formula)
        : remaining(formula.size()) {
    for (auto c : formula) {
      ClauseData cd;
      cd.literals = c;
      cd.orig_len = c.size();
      auto cpos = this->clauses.size();
      this->clauses.push_back(cd);
      for (auto l : c) this->add_literal(l, cpos);
    }
  }
};
/* formula ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::pure_literal][pure_literal]] */
void remove_satisfied(Formula& f, int d) {
  auto& clause = f.clauses[d];
  clause.n_t++;
  f.remaining--;
  if (VERBOSE) {
    std::cerr << "VERBOSE: Removing satisfied clause " << d << ", "
              << f.remaining << " Remaining" << std::endl;
  }
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
void pure_literal_assign(Formula& f, int l, Formula::LitData& data) {
  auto pos_size = data.pos_clauses.size();
  auto s = (pos_size == 0) ? -1 : 1;
  if (VERBOSE) std::cerr << "VERBOSE: Pure Literal - " << s*l << std::endl;
  auto lclauses = (s == 1) ? data.pos_clauses : data.neg_clauses;
  data.assn = (pos_size != 0) ? 1 : 0;
  if (VERBOSE) {
    std::cerr << "VERBOSE: Assigned " << l << " as "
              << (data.assn ? "TRUE" : "FALSE") << std::endl;
  }
  for (auto cidx : lclauses) remove_satisfied(f, cidx);
}
/* pure_literal ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::unit_propogate][unit_propogate]] */
void set_var(Formula& f, int l) {
  auto s = sign(l);
  auto pos = l*s;
  auto& lit = f.literals[pos];
  if (lit.assn != -1) throw std::runtime_error("literal already assigned");
  lit.assn = (s == 1) ? 1 : 0;
  if (VERBOSE) {
    std::cerr << "VERBOSE: Assigned " << l*s << " as "
              << (lit.assn ? "TRUE" : "FALSE") << std::endl;
  }
  auto sat_c = (lit.assn == 1) ? lit.pos_clauses : lit.neg_clauses;
  auto& unsat_c = (lit.assn == 0) ? lit.pos_clauses : lit.neg_clauses;
  for (auto cidx : sat_c) remove_satisfied(f, cidx);
  for (auto cidx : unsat_c) {
    auto& clause = f.clauses[cidx];
    clause.n_f++;
    if (VERBOSE) {
      std::cerr << "VERBOSE: Removing " << ((lit.assn == 0) ? pos : -pos)
                << " from clause " << cidx << std::endl;
    }
    clause.literals.erase(std::remove(clause.literals.begin(),
                                      clause.literals.end(),
                                      (lit.assn == 0) ? pos : -pos),
                          clause.literals.end());
  }
  unsat_c.clear();
}
void unit_propogate(Formula& f, Formula::ClauseData clause) {
  auto l = clause.literals[0];
  if (VERBOSE) std::cerr << "[" << l << "]" << std::endl;
  set_var(f, l);
}
/* unit_propogate ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::branch_helper][branch_helper]] */
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
  if (curr == 0) throw std::runtime_error("branching heuristic failed");
  return curr;
}
/* branch_helper ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::get_branching][get_branching]] */
auto dlis(Formula f, int l) {
  int wp = nclauses(f, -1, l);
  int wn = nclauses(f, -1, -l);
  return std::make_tuple(wp, wn, std::max(wp, wn));
}
auto dlcs(Formula f, int l) {
  int wp = nclauses(f, -1, l);
  int wn = nclauses(f, -1, -l);
  return std::make_tuple(wp, wn, wp + wn);
}
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
int get_branching_variable(Formula f, BranchRule rule) {
  switch (rule) {
    case BranchRule::dlis:
      return apply_rule(f, &dlis);
    case BranchRule::dlcs:
      return apply_rule(f, &dlcs);
    case BranchRule::jw:
      return apply_rule(f, &jw);
    case BranchRule::jw2:
      return apply_rule(f, &jw2);
    case BranchRule::dsj:
      return apply_rule(f, &dsj);
  }
  throw std::runtime_error("get_branching_variable didn't handle all cases");
}
/* get_branching ends here */
/* [[file:~/Documents/Projects/sat-solver/serial/serial.org::dpll][dpll]] */
std::tuple<bool, std::vector<int>, std::vector<int>> ret_val(Formula& f, bool sat) {
  if (!sat) return {false, {}, {}};
  std::vector<int> assnt;
  std::vector<int> assnf;
  for (auto l : f.literals) {
    if (l.second.assn == 1)
      assnt.push_back(l.first);
    if (l.second.assn == 0)
      assnf.push_back(l.first);
  }
  return {true, assnt, assnf};
}
std::tuple<bool, std::vector<int>, std::vector<int>> dpll(Formula& f, BranchRule rule) {
  for (auto& l : f.literals) {
    if (l.second.pure()) {
      pure_literal_assign(f, l.first, l.second);
    }
  }
  for (long unsigned int i = 0; i < f.clauses.size(); ++i) {
    auto c = f.clauses[i];
    if (c.sat()) continue;
    if (c.literals.size() == 0) {
      if (VERBOSE) std::cerr << "VERBOSE: Empty Clause - UNSAT" << std::endl;
      return ret_val(f, false);
    }
    if (c.unit()) {
      if (VERBOSE) std::cerr << "VERBOSE: Unit Propogate " << i << " ";
      unit_propogate(f, c);
    }
  }
  if (f.remaining == 0) {
    if (VERBOSE) std::cerr << "VERBOSE: No Remaining - SAT" << std::endl;
    return ret_val(f, true);
  }
  auto l = get_branching_variable(f, rule);
  if (VERBOSE) std::cerr << "VERBOSE: Branching on " << l << std::endl;
  Formula oldf(f);
  set_var(f, l);
  auto [res, ts, fs] = dpll(f, rule);
  if (res) return {res, ts, fs};
  if (VERBOSE) std::cerr << "VERBOSE: Backtracking" << std::endl;
  f = oldf;
  if (VERBOSE) std::cerr << "VERBOSE: Branching on " << -l << std::endl;
  set_var(f, -l);
  return dpll(f, rule);
}
/* dpll ends here */

void print_graph(Formula formula) {
  std::cerr << "VERBOSE: Literal graph is" << std::endl;
  std::cerr << "VERBOSE: LGRAPH digraph literal {" << std::endl;
  for (long unsigned int i = 0; i < formula.clauses.size(); ++i) {
    std::cerr << "VERBOSE: LGRAPH " << i << " [shape=box];" << std::endl;
  }
  for (auto l : formula.literals) {
    for (auto c : l.second.pos_clauses) {
      std::cerr << "VERBOSE: LGRAPH " << l.first << " -> " << c << ";" << std::endl;
    }
    for (auto c : l.second.neg_clauses) {
      std::cerr << "VERBOSE: LGRAPH " << -l.first << " -> " << c << ";" << std::endl;
    }
  }
  std::cerr << "VERBOSE: LGRAPH }" << std::endl;
}

void help(std::string name) {
  std::cout << "DPLL" << std::endl;
  std::cout << "Usage:\t" << name << " [OPTIONS]" << std::endl;
  std::cout << "\t-h\tShow this screen" << std::endl;
  std::cout << "\t-v\tVerbose Output" << std::endl;
  std::cout << "\t-r <RULE>\tSet branch rule (one of dlis|dlcs|jw|jw2|dsj)" << std::endl;
}

int main(int argc, char** argv) {
  int opt;
  std::string input = "";
  BranchRule rule = BranchRule::dlis;
  while ((opt = getopt(argc, argv, "hvr:")) != -1) {
    switch (opt) {
      case 'h': {
        help(argv[0]);
        return 0;
      }
      case 'v': {
        VERBOSE = true;
        break;
      }
      case 'r': {
        std::string r(optarg);
        if (r == "dlis") {
          rule = BranchRule::dlis;
        } else if (r == "dlcs") {
          rule = BranchRule::dlcs;
        } else if (r == "jw") {
          rule = BranchRule::jw;
        } else if (r == "jw2") {
          rule = BranchRule::jw2;
        } else if (r == "dsj") {
          rule = BranchRule::dsj;
        } else {
          help(argv[0]);
          return EXIT_FAILURE;
        }
        break;
      }
      default:
        help(argv[0]);
        return EXIT_FAILURE;
    }
  }
  
  Formula formula(read_input());
  if (VERBOSE) print_graph(formula);
  auto [sat, assnt, assnf] = dpll(formula, rule);
  std::cout << "Formula is: " << (sat ? "SAT" : "UNSAT") << std::endl;
  if (!sat) return 20;
  std::cout << "Variables assigned TRUE:" << std::endl;
  for (auto l : assnt) std::cout << l << " ";
  std::cout << std::endl;
  std::cout << "Variables assigned FALSE:" << std::endl;
  for (auto l : assnf) std::cout << l << " ";
  std::cout << std::endl;
  return 10;
}
/* Full Source:1 ends here */
