#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/bp_support.hpp>
#include <iostream>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <map>

template<class A, class B>
std::ostream& operator<<(std::ostream& os, const std::pair<A, B>& p) {
  return os << "(" << p.first << ", " << p.second << ")";
}

namespace sdsl {
  template<typename T, typename Idx>
  void construct_from_iv(Idx& idx, const T& iv) {
    int_vector<> v = iv;
    std::string tmp_file = ram_file_name(util::to_string(util::pid())+"_"+util::to_string(util::id()));
    store_to_file(v, tmp_file);
    construct(idx, tmp_file);
    ram_fs::remove(tmp_file);
  }

  class lz_forward_trie {
  public:
    struct node {
      int par;
      int par_edge;
      int offset;
      std::map<int,int> adj;
      node(int par=0, int par_edge=0, int offset=0)
        : par(par), par_edge(par_edge), offset(offset) {}
    };
    std::vector<node> nodes;

  private:
    void print_rec(std::ostream& out, int v, int indent=0) {
      out << v << " (offset=" << nodes[v].offset
              //<< ", par=" << nodes[v].par
              //<< ", par_edge=" << nodes[v].par_edge
              << ")" << std::endl;
      for (auto& it: nodes[v].adj) {
        for (int i = 0; i < indent + 1; ++i)
          out << "  ";
        out << (char)it.first << " -> ";
        print_rec(out, it.second, indent + 1);
      }
    }

  public:
    lz_forward_trie() : nodes(1) {}

    template<typename It>
    void lz_insert(It a, It b) {
      int cur = 0;
      int len = 0;
      int i = 0;
      for (auto c = a; c != b; ++c, ++i) {
        len++;
        auto it = nodes[cur].adj.find(*c);
        if (it != end(nodes[cur].adj)) {
          cur = it->second;
        } else {
          nodes.emplace_back(cur, *c, i - len + 1);
          nodes[cur].adj[*c] = nodes.size() - 1;
          cur = 0;
          len = 0;
        }
      }
    }

    void print(std::ostream& out) {
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

  class lz_rev_trie {
  public:
    struct edge {
      int forward_node, len;
    };
    struct node {
      int id;
      std::map<int,int> adj;
      std::vector<std::pair<edge, int>> edges;
      node(int id=0): id(id) {}
    };
    const lz_forward_trie& ft;
    std::vector<node> nodes;

  private:
    std::string get_edge(const edge& e) {
      std::string res;
      int i = e.forward_node;
      int len = e.len;
      while(len--) {
        res += ft.nodes[i].par_edge;
        i = ft.nodes[i].par;
      }
      return res;
    }

    void print_rec(std::ostream& out, int v, int indent=0) {
      out << v << " (" << nodes[v].id << ")" << std::endl;
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
          nodes[cur].edges.push_back(std::pair<edge,int>({id, len}, final_node));
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
          nodes.back().edges.push_back(std::pair<edge,int>({id, len}, final_node));
          nodes.back().edges.push_back(std::pair<edge,int>({j, l2}, nxt));
          nodes.back().adj[ft.nodes[id].par_edge] = 0;
          nodes.back().adj[ft.nodes[j].par_edge] = 1;

          nxt = nodes.size() - 1;
          e.len = common;
          return;
        }
        if (l2) { // complete remaining suffix is a prefix of the edge label
          nodes.emplace_back(orig_id);
          int final_node = nodes.size() - 1;

          nodes[final_node].edges.push_back(std::pair<edge,int>({j, l2}, nxt));
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
    lz_rev_trie(const lz_forward_trie& ft) : ft(ft), nodes(1) {
      for (int i = 0; i < ft.nodes.size(); ++i) {
        //print(cout);
        //cout << "============" << endl;
        insert_forward_leaf(i);
      }
    }

    void print(std::ostream& out) {
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

  class lz_index {
  public:
    using size_type = int;
    bit_vector ft_bp, rt_bp;
    bp_support_sada<> ft_bp_support, rt_bp_support;
    rank_support_v5<10,2> rt_bp_rank10;
    select_support_mcl<10,2> rt_bp_select10;
    int_vector<> ft_edge, ft_ids, ft_id_to_tree, ft_offset, rt_ids;
    wt_int<> succ_wt;
    int invalid_id;

    lz_index(const lz_forward_trie& ft, const lz_rev_trie& rt)
      : ft_bp(ft.nodes.size()*2, 0)
      , rt_bp(rt.nodes.size()*2, 0)
      , ft_edge(ft.nodes.size())
      , ft_ids(ft.nodes.size())
      , ft_id_to_tree(ft.nodes.size())
      , ft_offset(ft.nodes.size())
      , rt_ids(rt.nodes.size())
      , invalid_id(ft.nodes.size()*2)
    {
      int paren_idx, node_idx;
      paren_idx = node_idx = 0;
      //cout << "lztrie pre-order " << endl;
      ft.dfs(
        [&](int v, int par_edge) {
          //cout << v << " ";
          ft_bp[paren_idx] = 1;
          ft_edge[node_idx] = (par_edge >= 0) ? par_edge : 0;
          ft_id_to_tree[v] = paren_idx;
          ft_ids[node_idx] = v;
          ft_offset[node_idx] = ft.nodes[v].offset;
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
      int_vector<> succ_iv(rt.nodes.size());
      rt.dfs(
        [&](int v, const lz_rev_trie::edge& par_edge) {
          //cout << v << " ";
          rt_bp[paren_idx] = 1;
          int id = rt.nodes[v].id;
          rt_ids[node_idx] = (id != -1) ? ft_id_to_tree[id] : invalid_id;
          succ_iv[node_idx] =
              (id != -1 && id + 1 < ft_id_to_tree.size())
                ? ft_id_to_tree[id + 1]
                : invalid_id;
          paren_idx++;
          node_idx++;
        },
        [&](int v, const lz_rev_trie::edge& par_edge) {
          rt_bp[paren_idx++] = 0;
        });
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
      util::assign(rt_bp_support, bp_support_sada<>(&rt_bp));
      util::init_support(rt_bp_rank10, &rt_bp);
      util::init_support(rt_bp_select10, &rt_bp);
      construct_from_iv<>(succ_wt, succ_iv);
      util::bit_compress(ft_edge);
      util::bit_compress(ft_ids);
      util::bit_compress(ft_id_to_tree);
      util::bit_compress(ft_offset);
      util::bit_compress(rt_ids);
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
      if (b1 == invalid_id) {
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
        if (x == invalid_id)
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
        if (std::min(x, y) == -1)
          continue;
        //cout << "i y x " << i << " " << y << " " << x << endl;
        int a1 = rt_bp_support.rank(y) - 1,
            a2 = rt_bp_support.rank(rt_bp_support.find_close(y)) - 1,
            b1 = x,
            b2 = ft_bp_support.find_close(x);
        //cout << "a1 a2 b1 b2 = " << a1 << " " << a2 << " " << b1 << " " << b2 << endl;
        for (auto point : succ_wt.range_search_2d(a1, a2, b1, b2).second)
          report(trie_offset(point.second) - i);
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

    int serialize(std::ostream& out, structure_tree_node* v, std::string name) const {
      structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
      int written_bytes = 0;
#define SERIALIZE_FIELD(n) written_bytes += n.serialize(out, child, #n)
      SERIALIZE_FIELD(ft_bp);
      SERIALIZE_FIELD(ft_bp_support);
      SERIALIZE_FIELD(ft_edge);
      SERIALIZE_FIELD(ft_ids);
      SERIALIZE_FIELD(ft_id_to_tree);
      SERIALIZE_FIELD(ft_offset);
      SERIALIZE_FIELD(rt_bp);
      SERIALIZE_FIELD(rt_bp_support);
      SERIALIZE_FIELD(rt_bp_rank10);
      SERIALIZE_FIELD(rt_bp_select10);
      SERIALIZE_FIELD(rt_ids);
      SERIALIZE_FIELD(succ_wt);
      structure_tree::add_size(child, written_bytes);
      return written_bytes;
#undef SERIALIZE_FIELD
    }
  };
}
