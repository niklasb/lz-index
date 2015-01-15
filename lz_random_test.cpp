#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "lz_index.hpp"
using namespace std;

vector<string> test_strings = {
  "mississippihippi$",
};

int main() {
  // generate some random test strings
  for (int i = 0; i < 1000; ++i) {
    int n = rand() % 50;
    int a = rand() % 26 + 1;
    string s;
    for (int j = 0; j < n; ++j)
      s += 'a' + (rand() % a);
    s += '$';
    test_strings.push_back(s);
  }
  // test against brute force
  for (auto s: test_strings) {
    cout << "testing " << s << endl;
    sdsl::forward_trie ft;
    ft.lz_insert(begin(s), end(s));
    sdsl::rev_trie rt(ft);
    sdsl::succinct_lzindex lz(ft, rt);

    map<string, vector<int>> substrs;
    for (int i = 0; i < s.size(); ++i)
      for (int j = i; j < s.size(); ++j)
        substrs[s.substr(i, j - i + 1)].push_back(i);
    for (auto it: substrs) {
      sort(begin(it.second), end(it.second));
      vector<int> matches;
      lz.search_pattern(begin(it.first), end(it.first), [&](int i) {
        matches.push_back(i);
      });
      sort(begin(matches), end(matches));
      assert(matches == it.second);
    }
  }
}
