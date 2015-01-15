#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/bp_support.hpp>
#include <iostream>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <map>
using namespace std;
using namespace sdsl;

template<class A, class B>
ostream& operator<<(ostream& os, const std::pair<A, B>& p) {
  return os << "(" << p.first << ", " << p.second << ")";
}

template<typename T, typename Idx>
void construct_from_iv(Idx& idx, const T& iv) {
  int_vector<> v = iv;
  std::string tmp_file = ram_file_name(util::to_string(util::pid())+"_"+util::to_string(util::id()));
  store_to_file(v, tmp_file);
  construct(idx, tmp_file);
  ram_fs::remove(tmp_file);
}

class forward_trie {
public:
  struct node {
    int par;
    int par_edge;
    int offset;
    unordered_map<int,int> adj;
    node(int par=0, int par_edge=0, int offset=0)
      : par(par), par_edge(par_edge), offset(offset) {}
  };
  vector<node> nodes;

private:
  void print_rec(ostream& out, int v, int indent=0) {
    out << v << " (offset=" << nodes[v].offset
             //<< ", par=" << nodes[v].par
             //<< ", par_edge=" << nodes[v].par_edge
             << ")" << endl;
    for (auto& it: nodes[v].adj) {
      for (int i = 0; i < indent + 1; ++i)
        out << "  ";
      out << (char)it.first << " -> ";
      print_rec(out, it.second, indent + 1);
    }
  }

public:
  forward_trie() : nodes(1) {}

  template<typename S>
  void lz_insert(const S& s) {
    int cur = 0;
    int len = 0;
    for (int i = 0; i < s.size(); ++i) {
      len++;
      auto it = nodes[cur].adj.find(s[i]);
      if (it != end(nodes[cur].adj)) {
        cur = it->second;
      } else {
        nodes.emplace_back(cur, s[i], i - len + 1);
        nodes[cur].adj[s[i]] = nodes.size() - 1;
        cur = 0;
        len = 0;
      }
    }
  }

  void print(ostream& out) {
    print_rec(out, 0);
  }

  template <typename F, typename G>
  void dfs_(int v, F pre, G post) const {
    for (auto& it: nodes[v].adj) {
      pre(it.second, it.first);
      dfs_(it.second, pre, post);
      post(it.second, it.first);
    }
  }

  template <typename F, typename G>
  void dfs(F pre, G post) const {
    pre(0, -1);
    dfs_(0, pre, post);
    post(0, -1);
  }
};

class rev_trie {
public:
  struct edge {
    int forward_node, len;
  };
  struct node {
    int id;
    unordered_map<int,int> adj;
    vector<pair<edge, int>> edges;
    node(int id=0): id(id) {}
  };
  const forward_trie& ft;
  vector<node> nodes;

private:
  string get_edge(const edge& e) {
    string res;
    int i = e.forward_node;
    int len = e.len;
    while(len--) {
      res += ft.nodes[i].par_edge;
      i = ft.nodes[i].par;
    }
    return res;
  }

  void print_rec(ostream& out, int v, int indent=0) {
    out << v << " (" << nodes[v].id << ")" << endl;
    for (auto& e: nodes[v].edges) {
      for (int i = 0; i < indent + 1; ++i)
        out << "  ";
      out << get_edge(e.first) << " -> ";
      print_rec(out, e.second, indent + 1);
    }
  }

  void insert_forward_leaf(int id) {
    //cout<< "inserting " << id << endl;
    int len = 0;
    for (int i = id; i; i = ft.nodes[i].par) {
      //cout << "i = " << i << " par = " << ft.nodes[i].par << endl;
      len++;
    }
    int orig_id = id;
    //cout << "len = " << len << endl;
    int cur = 0;
    while (id) {
      //cout << "cur id len = "  << cur << " " << id << " " << len << endl;
      auto it = nodes[cur].adj.find(ft.nodes[id].par_edge);
      if (it == end(nodes[cur].adj)) {
        nodes.emplace_back(orig_id);
        int final_node = nodes.size() - 1;
        nodes[cur].edges.push_back(pair<edge,int>({id, len}, final_node));
        nodes[cur].adj[ft.nodes[id].par_edge] = nodes[cur].edges.size() - 1;
        return;
      }
      edge& e = nodes[cur].edges[it->second].first;
      int& nxt = nodes[cur].edges[it->second].second;
      int j = e.forward_node;
      int l2 = e.len;
      int common = 0;
      while (len && l2 && ft.nodes[id].par_edge == ft.nodes[j].par_edge) {
        id = ft.nodes[id].par;
        j = ft.nodes[j].par;
        len--;
        l2--;
        common++;
      }
      //cout<< "common = " << common << endl;
      if (!len && !l2) { // fill a node that was a split node before
        assert(nodes[nxt].id == -1);
        nodes[nxt].id = orig_id;
        return;
      }
      if (len && l2) { // we need to introduce a split node
        nodes.emplace_back(orig_id);
        int final_node = nodes.size() - 1;

        nodes.emplace_back(-1);
        nodes.back().edges.push_back(pair<edge,int>({id, len}, final_node));
        nodes.back().edges.push_back(pair<edge,int>({j, l2}, nxt));
        nodes.back().adj[ft.nodes[id].par_edge] = 0;
        nodes.back().adj[ft.nodes[j].par_edge] = 1;

        nxt = nodes.size() - 1;
        e.len = common;
        return;
      }
      if (l2) { // complete remaining suffix is a prefix of the edge label
        nodes.emplace_back(orig_id);
        int final_node = nodes.size() - 1;

        nodes[final_node].edges.push_back(pair<edge,int>({j, l2}, nxt));
        nodes[final_node].adj[ft.nodes[j].par_edge] = 0;
        nxt = final_node;
        e.len = common;
        return;
      }
      // edge label is prefix of the remaining suffix, recurse
      cur = nxt;
    }
  }

public:
  rev_trie(const forward_trie& ft) : ft(ft), nodes(1) {
    for (int i = 0; i < ft.nodes.size(); ++i) {
      //print(cout);
      //cout << "============" << endl;
      insert_forward_leaf(i);
    }
  }

  void print(ostream& out) {
    print_rec(out, 0);
  }

  template <typename F, typename G>
  void dfs_(int v, F pre, G post) const {
    for (auto& it: nodes[v].edges) {
      pre(it.second, it.first);
      dfs_(it.second, pre, post);
      post(it.second, it.first);
    }
  }

  template <typename F, typename G>
  void dfs(F pre, G post) const {
    pre(0, (edge){-1,-1});
    dfs_(0, pre, post);
    post(0, (edge){-1,-1});
  }
};

class succinct_lzindex {
public:
  bit_vector ft_bp, rt_bp;
  bp_support_sada<> ft_bp_support, rt_bp_support;
  rank_support_v5<10,2> ft_bp_rank10, rt_bp_rank10;
  select_support_mcl<10,2> ft_bp_select10, rt_bp_select10;
  int_vector<> ft_edge, ft_ids, ft_id_to_tree, ft_offset, rt_ids, succ_iv;
  wt_int<> succ_wt;
  int no_succ;

  succinct_lzindex(const forward_trie& ft, const rev_trie& rt)
    : ft_bp(ft.nodes.size()*2, 0)
    , rt_bp(rt.nodes.size()*2, 0)
    , ft_edge(ft.nodes.size())
    , ft_ids(ft.nodes.size())
    , ft_id_to_tree(ft.nodes.size())
    , ft_offset(ft.nodes.size())
    , rt_ids(rt.nodes.size())
    , succ_iv(rt.nodes.size() + 1)
    , no_succ(ft.nodes.size())
  {
    int paren_idx, node_idx;
    paren_idx = node_idx = 0;
    //cout << "lztrie pre-order " << endl;
    vector<int> pos_id(ft.nodes.size());
    ft.dfs(
      [&](int v, int par_edge) {
        //cout << v << " ";
        ft_bp[paren_idx] = 1;
        ft_edge[node_idx] = par_edge;
        ft_id_to_tree[v] = paren_idx;
        ft_ids[node_idx] = v;
        ft_offset[node_idx] = ft.nodes[v].offset;
        pos_id[v] = node_idx;
        paren_idx++;
        node_idx++;
      },
      [&](int v, int par_edge) {
        ft_bp[paren_idx++] = 0;
      });
    //cout << endl;
    //cout << "trie seq" << endl;
    //for (int i = 0; i < ft_bp.size(); ++i)
      //cout << ft_bp[i] << " ";
    //cout << endl;
    paren_idx = node_idx = 0;
    //cout << ft_id_to_tree.size() << endl;
    //cout << "revtrie pre-order " << endl;
    rt.dfs(
      [&](int v, const rev_trie::edge& par_edge) {
        //cout << v << " ";
        rt_bp[paren_idx] = 1;
        int id = rt.nodes[v].id;
        rt_ids[node_idx] = (id != -1) ? ft_id_to_tree[id] : -1;
        if (id != -1 && id + 1 < pos_id.size()) {
          succ_iv[node_idx] = pos_id[id + 1];
          assert(succ_iv[node_idx] != no_succ);
        } else {
          succ_iv[node_idx] = no_succ;
        }
        paren_idx++;
        node_idx++;
      },
      [&](int v, const rev_trie::edge& par_edge) {
        rt_bp[paren_idx++] = 0;
      });
    succ_iv[succ_iv.size() - 1] = no_succ;
    assert(succ_iv.size() > no_succ);
    //cout << endl;
    //cout << "revtrie seq" << endl;
    //for (int i = 0; i < rt_bp.size(); ++i)
      //cout << rt_bp[i] << " ";
    //cout << endl;
    //cout << "succ = " << endl;
    //for (int i = 0; i < succ_iv.size(); ++i)
      //cout << succ_iv[i] << " ";
    //cout << endl;
    util::assign(ft_bp_support, bp_support_sada<>(&ft_bp));
    util::init_support(ft_bp_rank10, &ft_bp);
    util::init_support(ft_bp_select10, &ft_bp);
    util::assign(rt_bp_support, bp_support_sada<>(&rt_bp));
    util::init_support(rt_bp_rank10, &rt_bp);
    util::init_support(rt_bp_select10, &rt_bp);
    construct_from_iv<>(succ_wt, succ_iv);
  }

  int trie_edge(int x) {
    return ft_edge[ft_bp_support.rank(x) - 1];
  }

  int trie_offset(int x) {
    return ft_offset[ft_bp_support.rank(x) - 1];
  }

  int trie_id(int x) {
    return ft_ids[ft_bp_support.rank(x) - 1];
  }

  int trie_child(int x, int c) {
    int child = x + 1;
    while (ft_bp[child]) {
      if (trie_edge(child) == c)
        return child;
      child = ft_bp_support.find_close(child) + 1;
    }
    return -1;
  }

  int trie_depth(int x) {
    int cnt = 0;
    while (x) {
      cnt++;
      x = ft_bp_support.enclose(x);
    }
    return cnt;
  }

  template<typename It>
  int trie_search(int x, It a, It b) {
    for (auto it = a; it != b; ++it) {
      //cout << "trie_child(" << x << ", " << *it << ") = " << trie_child(x, *it) << endl;
      x = trie_child(x, *it);
      if (x == -1) return -1;
    }
    return x;
  }

  template<typename It>
  bool trie_endswith(int x, It a, It b) {
    if (a == b) return true;
    auto it = b;
    while (it != a) {
      --it;
      if (!x || *it != trie_edge(x))
        return false;
      x = ft_bp_support.enclose(x);
    }
    return true;
  }

  template<typename It>
  bool trie_startswith(int x, It a, It b) {
    int sz = trie_depth(x);
    //cout << "sz = " << sz << endl;
    if (sz < b - a) return false;
    if (a == b) return true;
    //cout << "skip = " << (sz-(b-a)) << endl;
    for (int i = 0; i < sz - (b - a); ++i)
      x = ft_bp_support.enclose(x);
    //cout << "x = " << x << endl;
    return trie_endswith(x, a, b);
  }

  int revtrie_to_trie(int y) {
    return rt_ids[rt_bp_support.rank(y) - 1];
  }

  // i == 1-based
  int revtrie_select_leaf(int i) {
    return rt_bp_select10.select(i) - 1;
  }

  int revtrie_depth(int y) {
    //cout << "y = " << y << endl;
    int rank_first_leaf = rt_bp_rank10(y) + 1;
    //cout << "rank_first_leaf = " << rank_first_leaf << endl;
    //cout << "b0 revtrie = " << revtrie_select_leaf(rank_first_leaf) << endl;
    int b0 = revtrie_to_trie(revtrie_select_leaf(rank_first_leaf));
    int b1 = revtrie_to_trie(y);
    if (b1 == -1) {
      // there has to be a second leaf
      //cout << "b1 revtrie = " << revtrie_select_leaf(rank_first_leaf + 1) << endl;
      b1 = revtrie_to_trie(revtrie_select_leaf(rank_first_leaf + 1));
    }
    //cout << "b0 b1 = " << b0 << " " << b1 << endl;
    int cnt = 0;
    while (b0 && b1 && trie_edge(b0) == trie_edge(b1)) {
      cnt++;
      b0 = ft_bp_support.enclose(b0);
      b1 = ft_bp_support.enclose(b1);
    }
    return cnt;
  }

  // return subtree that has s as a prefix.
  // returns -1 for empty subtree
  template<typename It>
  int revtrie_search(It a, It b) {
    int s_offset = 0;
    int cur = 0;
    int n = b - a;
    while (s_offset < n) {
      //cout << "cur s_offset s = " << cur << " " << s_offset << " " << s<< endl;
      bool found = 0;
      int child = cur + 1;
      while (rt_bp[child]) {
        //cout << "  child = " << child << endl;
        // get forward trie node for first leaf in subtree
        int rank_first_leaf = rt_bp_rank10(child) + 1;
        int x = revtrie_to_trie(revtrie_select_leaf(rank_first_leaf));
        //cout << "  x = " << x << endl;
        // walk to split point
        for (int i = 0; i < s_offset; ++i)
          x = ft_bp_support.enclose(x);
        //cout << "  x = " << x << endl;
        // check next character
        if (trie_edge(x) == *(b - s_offset - 1)) {
          int len_diff = revtrie_depth(child) - revtrie_depth(cur);
          //cout << "  len_diff = " << len_diff << endl;
          for (int i = 0; i < len_diff; ++i) {
            //cout << "    " << (char)trie_edge(x) << endl;
            if (s_offset + i < n && trie_edge(x) != *(b - s_offset - i - 1)) {
              //cout << "    mismatch at " << i << endl;
              return -1;
            }
            x = ft_bp_support.enclose(x);
          }
          s_offset += len_diff;
          cur = child;
          found = 1;
          break;
        }
        child = rt_bp_support.find_close(child) + 1;
      }
      if (!found)
        return -1;
    }
    return cur;
  }

  template<typename T, typename F>
  void dfs(T bp, int x, F f) {
    int dep = 0;
    do {
      if (bp[x]) f(x, dep);
      dep += (bp[x] ? 1 : -1);
      ++x;
    } while (dep);
  }

  template<typename It, typename F>
  void case1(It a, It b, F report) {
    int subtree = revtrie_search(a, b);
    int n = b - a;
    if (subtree == -1) return;
    dfs(rt_bp, subtree, [&](int y, int _) {
      int x = revtrie_to_trie(y);
      if (x == -1)
        return;
      dfs(ft_bp, x, [&](int x, int depth_x) {
        report(trie_offset(x) + trie_depth(x) - n - depth_x);
      });
    });
  }

  template<typename It, typename F>
  void case2(It a, It b, F report) {
    int n = b - a;
    for (int i = 1; i < n; ++i) {
      int y = revtrie_search(a, a + i);
      int x = trie_search(0, a + i, b);
      if (min(x, y) == -1)
        continue;
      //cout << "i y x " << i << " " << y << " " << x << endl;
      int a1 = rt_bp_support.rank(y) - 1,
          a2 = rt_bp_support.rank(rt_bp_support.find_close(y)) - 1,
          b1 = ft_bp_support.rank(x) - 1,
          b2 = ft_bp_support.rank(ft_bp_support.find_close(x)) - 1;
      //cout << "a1 a2 b1 b2 = " << a1 << " " << a2 << " " << b1 << " " << b2 << endl;
      auto res = succ_wt.range_search_2d(a1, a2, b1, b2);
      for (auto point : res.second) {
        //cout << point << endl;
        report(ft_offset[point.second] - i);
      }
    }
  }

  template<typename It, typename F>
  void case3(It a, It b, F report) {
    int n = b - a;
    for (int l = 0; l < n; ++l) {
      int mid = 0;
      for (int r = l; r < n; ++r) {
        // we are matching the substring pattern[l..r]
        mid = trie_child(mid, *(a + r));
        if (mid == -1)
          break;
        //cout << "l r mid = " << l << " " << r << " " << mid << endl;

        int id = trie_id(mid);
        int covered_blocks = 1 + (l > 0);
        int prev_x, prev_sz, i;
        //cout << "id covered_blocks = " << id << " " << covered_blocks << endl;

        // to check:
        // (1) the block before mid needs to be partially (not fully!) covered from
        //     the end
        // (2) at least 3 blocks need to be covered

        if (id == 1 && l > 0)
          // there is no previous block but pattern continues towards the left
          goto fail;
        if (id > 1) {
          prev_x = ft_id_to_tree[id - 1];
          prev_sz = trie_depth(prev_x);
          if (l >= prev_sz)
            // previous block is fully covered
            goto fail;
          if (!trie_endswith(prev_x, a, a + l))
            // previous block does not end with pattern start
            goto fail;
        }

        // check to the right of mid
        ++id;
        i = r + 1;
        while (i < n) {
          covered_blocks++;
          if (id >= ft_id_to_tree.size())
            goto fail;
          int x = ft_id_to_tree[id];
          int sz = trie_depth(x);
          //cout << "id i x sz = " << id << " " << i << " " << x << " " << sz << endl;
          if (sz <= n - i) {
            // pattern covers whole block
            if (!trie_endswith(x, a + i, a + i + sz)) {
              //cout << "id " << id << " " << x << " does not end with " << i << ".." << (i+sz) << endl;
              goto fail;
            }
            id++;
            i += sz;
          } else {
            // pattern covers beginning of block
            if (!trie_startswith(x, a + i, b)) {
              //cout << "id " << id << " " << x << " does not start with " << i << "..-1" << endl;
              goto fail;
            }
            break;
          }
        }

        if (covered_blocks < 3)
          goto fail;

        // we found a match!
        report(trie_offset(mid) - l);
fail:;
      }
    }
  }

  template<typename It, typename F>
  void search_pattern(It a, It b, F report) {
    case1(a, b, report);
    case2(a, b, report);
    case3(a, b, report);
  }
};

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
    forward_trie ft;
    ft.lz_insert(s);
    rev_trie rt(ft);
    succinct_lzindex lz(ft, rt);

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
