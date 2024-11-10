#include "seq_io.hpp"

#include <cassert>
#include <algorithm>
#include <cassert>

namespace fsa {

bool FastaReader::Next(Item &item) {
    item.id = Tell();
    item.quality = "";
    return GetHead(item.head, item.sub_head) && GetSeq(item.seq);
}

bool FastaReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = in_.QueryNoEmptyLine();
    in_.ConsumeStrippedLine();
    //printf("xxx：%s\n", line.c_str()); fflush(stdout);
    if (!line.empty() && line[0] == '>') {
        std::string::size_type s = 1;
        while (s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s+1, line.size());
        while (e < line.size() && !::isspace(line[e])) e++;

        if (e > s) {
            head = line.substr(s, e-s);
            while (e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool FastaReader::GetSeq(std::string &seq) {
    seq = "";
    std::string line = in_.QueryNoEmptyLine();

    while (!line.empty() && line[0] != '>') {
        seq += line;
        in_.ConsumeStrippedLine();
        line = in_.QueryNoEmptyLine();
    }

    return true;// !seq.empty(); allow empty sequence

}



bool FastqReader::Next(Item &item) {
    assert(IsValid());

    item.id = Tell();
    return GetHead(item.head, item.sub_head) && GetSeq(item.seq) &&
           GetHead1() && GetQuality(item.quality);
}

bool FastqReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = in_.GetStrippedLine();

    if (!line.empty() && line[0] == '@') {
        std::string::size_type s = 1;
        while (s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s+1, line.size());
        while (e < line.size() && !::isspace(line[e])) e++;

        if (e > s) {
            head = line.substr(s, e-s);
            while (e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

//////////////////
bool FastaInBlock::Next(Item& item) {
    item.quality[0] = item.quality[1] = nullptr;
    return GetHead(item) && GetSeq(item);

}
bool FastaInBlock::GetHead(Item& item) {
    index_ = std::find(index_, end_, '>');
    if (index_ < end_) {
        auto s = index_ + 1;
        while (s < end_ && ::isspace(*s)) s++;
        auto e = s;
        while (e < end_ && ::isspace(*e)) e++;
        if (e > s) {
            item.head[0] = s;
            item.head[1] = e;

            auto ss = e;
            while (ss < end_ && ::isspace(*ss)) ss++;
            auto se = ss;
            while (se < end_ && *se != '\n') se++;
            item.sub_head[0] = ss;
            item.sub_head[1] = se;
            return true;
        } else {
            return false;

        }

    } else {
        return false;
    }
}

bool FastaInBlock::GetSeq(Item& item) {
    
    return true;
}


} // namespace fsa {
    