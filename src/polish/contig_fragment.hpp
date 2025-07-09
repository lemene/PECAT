namespace fsa {
class ContigFragment {
public:
    ContigFragment(Seq::Id id, const DnaSeq& seq, size_t start, size_t end)
        : id_(id), seq_(seq), start_(start), end_(end) {}       

}

}