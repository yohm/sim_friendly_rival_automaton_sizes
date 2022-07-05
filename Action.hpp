#ifndef ACTION_HPP
#define ACTION_HPP


enum Action {
  C,
  D
};

char A2C(Action act);
Action C2A(char c);

std::ostream &operator<<(std::ostream &os, const Action &act);
#endif // ACTION_HPP
