#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <iterator>
#include "lz_index.hpp"

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " Creates an LZ index for a byte file and visualizes the memory utilization during construction." << endl;
        return 1;
    }

    memory_monitor::start();

    auto start = timer::now();
    ifstream f(argv[1]);
    f >> std::noskipws;
    sdsl::forward_trie ft;
    ft.lz_insert(istream_iterator<char>(f), istream_iterator<char>());
    sdsl::rev_trie rt(ft);
    sdsl::succinct_lzindex lz(ft, rt);
    auto stop = timer::now();
    cout << "construction time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;

    memory_monitor::stop();

    std::cout << "peak usage = " << memory_monitor::peak() / (1024*1024) << " MB" << std::endl;

    std::ofstream cstofs("lz-construction.html");
    cout << "writing memory usage visualization to lz-construction.html\n";
    memory_monitor::write_memory_log<HTML_FORMAT>(cstofs);
    cstofs.close();
}
