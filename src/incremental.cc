#include "incremental.h"
#include "service.h"
#include "discovery.h"
#include "unifier.h"
#include "genotyper.h"
#include "cli_utils.h"
#include "RocksKeyValue.h"
#include "vcf.h"
#include "ctpl_stl.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cassert>
#include <cstring>
#include <atomic>
#include <endian.h>

using namespace std;

namespace GLnexus {

// ============================================================================
// Column family names
// ============================================================================

const char* DISCOVERED_CF = "discovered";
const char* UNIFIED_CF = "unified";
const char* GENOTYPE_CACHE_CF = "genotype_cache";
const char* INCREMENTAL_META_CF = "incremental_meta";

// ============================================================================
// Binary serialization helpers
// ============================================================================

class BinaryWriter {
    string buf_;
public:
    void write_int32(int32_t v) { buf_.append(reinterpret_cast<const char*>(&v), sizeof(v)); }
    void write_uint32(uint32_t v) { buf_.append(reinterpret_cast<const char*>(&v), sizeof(v)); }
    void write_float(float v) { buf_.append(reinterpret_cast<const char*>(&v), sizeof(v)); }
    void write_uint8(uint8_t v) { buf_.append(reinterpret_cast<const char*>(&v), 1); }
    void write_string(const string& s) {
        write_uint32(s.size());
        buf_.append(s);
    }
    const string& data() const { return buf_; }
};

class BinaryReader {
    const char* data_;
    size_t size_;
    size_t pos_ = 0;
public:
    BinaryReader(const char* data, size_t size) : data_(data), size_(size) {}

    Status read_int32(int32_t& v) {
        if (pos_ + sizeof(int32_t) > size_) return Status::Invalid("BinaryReader: read past end");
        memcpy(&v, data_ + pos_, sizeof(int32_t));
        pos_ += sizeof(int32_t);
        return Status::OK();
    }
    Status read_uint32(uint32_t& v) {
        if (pos_ + sizeof(uint32_t) > size_) return Status::Invalid("BinaryReader: read past end");
        memcpy(&v, data_ + pos_, sizeof(uint32_t));
        pos_ += sizeof(uint32_t);
        return Status::OK();
    }
    Status read_float(float& v) {
        if (pos_ + sizeof(float) > size_) return Status::Invalid("BinaryReader: read past end");
        memcpy(&v, data_ + pos_, sizeof(float));
        pos_ += sizeof(float);
        return Status::OK();
    }
    Status read_uint8(uint8_t& v) {
        if (pos_ + 1 > size_) return Status::Invalid("BinaryReader: read past end");
        v = static_cast<uint8_t>(data_[pos_]);
        pos_ += 1;
        return Status::OK();
    }
    Status read_string(string& s) {
        uint32_t len;
        Status st;
        st = read_uint32(len);
        if (st.bad()) return st;
        if (pos_ + len > size_) return Status::Invalid("BinaryReader: string read past end");
        s.assign(data_ + pos_, len);
        pos_ += len;
        return Status::OK();
    }
    bool at_end() const { return pos_ >= size_; }
};

// ============================================================================
// discovered_alleles serialization
// ============================================================================

// Binary format version
static const uint32_t DISCOVERED_ALLELES_VERSION = 1;

Status serialize_discovered_alleles(const discovered_alleles& dsals, string& out) {
    BinaryWriter w;
    w.write_uint32(DISCOVERED_ALLELES_VERSION);
    w.write_uint32(dsals.size());

    for (const auto& p : dsals) {
        const allele& al = p.first;
        const discovered_allele_info& info = p.second;

        // allele: range + dna
        w.write_int32(al.pos.rid);
        w.write_int32(al.pos.beg);
        w.write_int32(al.pos.end);
        w.write_string(al.dna);

        // info
        w.write_uint8(info.is_ref ? 1 : 0);
        w.write_uint8(info.all_filtered ? 1 : 0);

        // topAQ: 10 int values
        for (unsigned i = 0; i < top_AQ::COUNT; i++) {
            w.write_int32(info.topAQ.V[i]);
        }

        // zGQ: 10 x 2 unsigned values
        for (unsigned i = 0; i < zygosity_by_GQ::GQ_BANDS; i++) {
            for (unsigned j = 0; j < zygosity_by_GQ::PLOIDY; j++) {
                w.write_uint32(info.zGQ.M[i][j]);
            }
        }

        // in_target range
        w.write_int32(info.in_target.rid);
        w.write_int32(info.in_target.beg);
        w.write_int32(info.in_target.end);
    }

    out = w.data();
    return Status::OK();
}

Status deserialize_discovered_alleles(const char* data, size_t size, discovered_alleles& ans) {
    BinaryReader r(data, size);
    Status s;

    uint32_t version;
    S(r.read_uint32(version));
    if (version != DISCOVERED_ALLELES_VERSION) {
        return Status::Invalid("discovered_alleles cache version mismatch");
    }

    uint32_t count;
    S(r.read_uint32(count));

    ans.clear();
    for (uint32_t i = 0; i < count; i++) {
        int32_t rid, beg, end;
        S(r.read_int32(rid));
        S(r.read_int32(beg));
        S(r.read_int32(end));

        string dna;
        S(r.read_string(dna));

        allele al(range(rid, beg, end), dna);

        discovered_allele_info info;
        uint8_t b;
        S(r.read_uint8(b)); info.is_ref = (b != 0);
        S(r.read_uint8(b)); info.all_filtered = (b != 0);

        for (unsigned j = 0; j < top_AQ::COUNT; j++) {
            S(r.read_int32(info.topAQ.V[j]));
        }
        for (unsigned j = 0; j < zygosity_by_GQ::GQ_BANDS; j++) {
            for (unsigned k = 0; k < zygosity_by_GQ::PLOIDY; k++) {
                S(r.read_uint32(info.zGQ.M[j][k]));
            }
        }

        int32_t it_rid, it_beg, it_end;
        S(r.read_int32(it_rid));
        S(r.read_int32(it_beg));
        S(r.read_int32(it_end));
        info.in_target = range(it_rid, it_beg, it_end);

        ans[al] = info;
    }
    return Status::OK();
}

Status deserialize_discovered_alleles(const string& data, discovered_alleles& ans) {
    return deserialize_discovered_alleles(data.data(), data.size(), ans);
}

// ============================================================================
// unified_site serialization
// ============================================================================

static const uint32_t UNIFIED_SITES_VERSION = 1;

Status serialize_unified_sites(const vector<unified_site>& sites, string& out) {
    BinaryWriter w;
    w.write_uint32(UNIFIED_SITES_VERSION);
    w.write_uint32(sites.size());

    for (const auto& us : sites) {
        // pos
        w.write_int32(us.pos.rid);
        w.write_int32(us.pos.beg);
        w.write_int32(us.pos.end);

        // in_target
        w.write_int32(us.in_target.rid);
        w.write_int32(us.in_target.beg);
        w.write_int32(us.in_target.end);

        // alleles
        w.write_uint32(us.alleles.size());
        for (const auto& ua : us.alleles) {
            w.write_string(ua.dna);
            w.write_int32(ua.normalized.pos.rid);
            w.write_int32(ua.normalized.pos.beg);
            w.write_int32(ua.normalized.pos.end);
            w.write_string(ua.normalized.dna);
            w.write_int32(ua.quality);
            w.write_float(ua.frequency);
        }

        // unification map
        w.write_uint32(us.unification.size());
        for (const auto& p : us.unification) {
            w.write_int32(p.first.pos.rid);
            w.write_int32(p.first.pos.beg);
            w.write_int32(p.first.pos.end);
            w.write_string(p.first.dna);
            w.write_int32(p.second);
        }

        w.write_float(us.lost_allele_frequency);
        w.write_int32(us.qual);
        w.write_uint8(us.monoallelic ? 1 : 0);
    }

    out = w.data();
    return Status::OK();
}

Status deserialize_unified_sites(const char* data, size_t size, vector<unified_site>& ans) {
    BinaryReader r(data, size);
    Status s;

    uint32_t version;
    S(r.read_uint32(version));
    if (version != UNIFIED_SITES_VERSION) {
        return Status::Invalid("unified_sites cache version mismatch");
    }

    uint32_t num_sites;
    S(r.read_uint32(num_sites));

    ans.clear();
    ans.reserve(num_sites);

    for (uint32_t i = 0; i < num_sites; i++) {
        int32_t rid, beg, end;
        S(r.read_int32(rid));
        S(r.read_int32(beg));
        S(r.read_int32(end));
        unified_site us(range(rid, beg, end));

        S(r.read_int32(us.in_target.rid));
        S(r.read_int32(us.in_target.beg));
        S(r.read_int32(us.in_target.end));

        uint32_t num_alleles;
        S(r.read_uint32(num_alleles));
        for (uint32_t j = 0; j < num_alleles; j++) {
            string dna;
            S(r.read_string(dna));
            unified_allele ua(us.pos, dna);

            S(r.read_int32(ua.normalized.pos.rid));
            S(r.read_int32(ua.normalized.pos.beg));
            S(r.read_int32(ua.normalized.pos.end));
            S(r.read_string(ua.normalized.dna));
            S(r.read_int32(ua.quality));
            S(r.read_float(ua.frequency));

            us.alleles.push_back(ua);
        }

        uint32_t num_unification;
        S(r.read_uint32(num_unification));
        for (uint32_t j = 0; j < num_unification; j++) {
            int32_t u_rid, u_beg, u_end;
            S(r.read_int32(u_rid));
            S(r.read_int32(u_beg));
            S(r.read_int32(u_end));
            string u_dna;
            S(r.read_string(u_dna));
            int32_t to_idx;
            S(r.read_int32(to_idx));
            us.unification[allele(range(u_rid, u_beg, u_end), u_dna)] = to_idx;
        }

        S(r.read_float(us.lost_allele_frequency));
        S(r.read_int32(us.qual));
        uint8_t mono;
        S(r.read_uint8(mono));
        us.monoallelic = (mono != 0);

        ans.push_back(move(us));
    }
    return Status::OK();
}

Status deserialize_unified_sites(const string& data, vector<unified_site>& ans) {
    return deserialize_unified_sites(data.data(), data.size(), ans);
}

// ============================================================================
// Incremental collection helpers
// ============================================================================

string contig_key(int rid) {
    // 3-byte big-endian encoding (same as BCFBucketRange prefix for rid)
    uint64_t rid_be = htobe64(rid);
    char buf[3];
    memcpy(buf, ((char*)&rid_be) + 5, 3);
    return string(buf, 3);
}

Status ensure_incremental_collections(KeyValue::DB* db) {
    Status s;
    KeyValue::CollectionHandle coll;

    // Try to get each collection. If not found, create it.
    const char* incremental_collections[] = { DISCOVERED_CF, UNIFIED_CF, GENOTYPE_CACHE_CF, INCREMENTAL_META_CF };
    for (const auto& name : incremental_collections) {
        s = db->collection(name, coll);
        if (s == StatusCode::NOT_FOUND) {
            s = db->create_collection(name);
            if (s.bad() && s != StatusCode::EXISTS) {
                return s;
            }
        } else if (s.bad()) {
            return s;
        }
    }
    return Status::OK();
}

Status put_incremental_meta(KeyValue::DB* db, const string& key, const string& value) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(INCREMENTAL_META_CF, coll));
    return db->put(coll, key, value);
}

Status get_incremental_meta(KeyValue::DB* db, const string& key, string& value) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(INCREMENTAL_META_CF, coll));
    return db->get(coll, key, value);
}

Status put_cached_discovered_alleles(KeyValue::DB* db, int rid, const discovered_alleles& dsals) {
    Status s;
    string serialized;
    S(serialize_discovered_alleles(dsals, serialized));

    KeyValue::CollectionHandle coll;
    S(db->collection(DISCOVERED_CF, coll));
    return db->put(coll, contig_key(rid), serialized);
}

Status get_cached_discovered_alleles(KeyValue::DB* db, int rid, discovered_alleles& dsals) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(DISCOVERED_CF, coll));

    string serialized;
    S(db->get(coll, contig_key(rid), serialized));

    return deserialize_discovered_alleles(serialized, dsals);
}

Status put_cached_unified_sites(KeyValue::DB* db, int rid, const vector<unified_site>& sites) {
    Status s;
    string serialized;
    S(serialize_unified_sites(sites, serialized));

    KeyValue::CollectionHandle coll;
    S(db->collection(UNIFIED_CF, coll));
    return db->put(coll, contig_key(rid), serialized);
}

Status get_cached_unified_sites(KeyValue::DB* db, int rid, vector<unified_site>& sites) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(UNIFIED_CF, coll));

    string serialized;
    S(db->get(coll, contig_key(rid), serialized));

    return deserialize_unified_sites(serialized, sites);
}

// ============================================================================
// Phase 2: Genotype cache (per site)
// ============================================================================

string site_key(int rid, int beg, int end) {
    // 11-byte key: 3 bytes rid (BE) + 4 bytes beg (BE) + 4 bytes end (BE)
    uint64_t rid_be = htobe64(rid);
    uint32_t beg_be = htobe32(beg);
    uint32_t end_be = htobe32(end);
    char buf[11];
    memcpy(buf, ((char*)&rid_be) + 5, 3);
    memcpy(buf + 3, &beg_be, 4);
    memcpy(buf + 7, &end_be, 4);
    return string(buf, 11);
}

Status put_cached_genotype(KeyValue::DB* db, const range& pos, const string& bcf_data) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(GENOTYPE_CACHE_CF, coll));
    return db->put(coll, site_key(pos.rid, pos.beg, pos.end), bcf_data);
}

Status get_cached_genotype(KeyValue::DB* db, const range& pos, string& bcf_data) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(GENOTYPE_CACHE_CF, coll));
    return db->get(coll, site_key(pos.rid, pos.beg, pos.end), bcf_data);
}

Status clear_genotype_cache(KeyValue::DB* db) {
    Status s;
    KeyValue::CollectionHandle coll;
    S(db->collection(GENOTYPE_CACHE_CF, coll));

    // Iterate and collect all keys, then delete by overwriting with empty
    // (RocksDB doesn't expose delete through our KeyValue interface,
    //  so we overwrite with empty values as a soft clear)
    unique_ptr<KeyValue::Iterator> it;
    S(db->iterator(coll, "", it));
    unique_ptr<KeyValue::WriteBatch> wb;
    S(db->begin_writes(wb));
    size_t count = 0;
    while (it->valid()) {
        S(wb->put(coll, it->key().str(), ""));
        count++;
        S(it->next());
    }
    if (count > 0) {
        S(wb->commit());
    }
    return Status::OK();
}

// ============================================================================
// Phase 3: Validation — compare discovered alleles from incremental vs full
// ============================================================================

Status validate_bcf_genotypes(const string& file_a, const string& file_b, string& report) {
    // Open both BCF files and compare GT fields record by record
    unique_ptr<vcfFile, void(*)(vcfFile*)> va(bcf_open(file_a.c_str(), "r"),
                                               [](vcfFile* f) { if(f) bcf_close(f); });
    unique_ptr<vcfFile, void(*)(vcfFile*)> vb(bcf_open(file_b.c_str(), "r"),
                                               [](vcfFile* f) { if(f) bcf_close(f); });
    if (!va || !vb) {
        return Status::IOError("validate_bcf_genotypes: failed to open BCF files");
    }

    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> ha(bcf_hdr_read(va.get()), &bcf_hdr_destroy);
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hb(bcf_hdr_read(vb.get()), &bcf_hdr_destroy);
    if (!ha || !hb) {
        return Status::IOError("validate_bcf_genotypes: failed to read headers");
    }

    // Compare sample counts
    int na = bcf_hdr_nsamples(ha.get());
    int nb = bcf_hdr_nsamples(hb.get());
    if (na != nb) {
        report = "Sample count mismatch: " + to_string(na) + " vs " + to_string(nb);
        return Status::Invalid("validation failed", report);
    }

    unique_ptr<bcf1_t, void(*)(bcf1_t*)> ra(bcf_init(), &bcf_destroy);
    unique_ptr<bcf1_t, void(*)(bcf1_t*)> rb(bcf_init(), &bcf_destroy);

    size_t sites_compared = 0, sites_matched = 0, sites_mismatched = 0;
    ostringstream mismatch_details;

    while (true) {
        int ret_a = bcf_read(va.get(), ha.get(), ra.get());
        int ret_b = bcf_read(vb.get(), hb.get(), rb.get());

        if (ret_a < 0 && ret_b < 0) break; // both done
        if (ret_a < 0 || ret_b < 0) {
            report = "Different number of records";
            return Status::Invalid("validation failed", report);
        }

        bcf_unpack(ra.get(), BCF_UN_ALL);
        bcf_unpack(rb.get(), BCF_UN_ALL);
        sites_compared++;

        // Compare position
        if (ra->rid != rb->rid || ra->pos != rb->pos) {
            mismatch_details << "Position mismatch at record " << sites_compared
                             << ": " << ra->rid << ":" << ra->pos
                             << " vs " << rb->rid << ":" << rb->pos << "\n";
            sites_mismatched++;
            continue;
        }

        // Compare GT fields
        int32_t *gt_a = nullptr, *gt_b = nullptr;
        int ngt_a = 0, ngt_b = 0;
        bcf_get_genotypes(ha.get(), ra.get(), &gt_a, &ngt_a);
        bcf_get_genotypes(hb.get(), rb.get(), &gt_b, &ngt_b);

        bool gt_match = (ngt_a == ngt_b);
        if (gt_match && ngt_a > 0) {
            gt_match = (memcmp(gt_a, gt_b, ngt_a * sizeof(int32_t)) == 0);
        }

        if (gt_match) {
            sites_matched++;
        } else {
            sites_mismatched++;
            if (sites_mismatched <= 10) {
                mismatch_details << "GT mismatch at " << ra->rid << ":" << ra->pos + 1 << "\n";
            }
        }

        free(gt_a);
        free(gt_b);
    }

    ostringstream summary;
    summary << "Compared " << sites_compared << " sites: "
            << sites_matched << " matched, " << sites_mismatched << " mismatched";
    if (sites_mismatched > 0) {
        summary << "\n" << mismatch_details.str();
    }
    report = summary.str();

    if (sites_mismatched > 0) {
        return Status::Invalid("validation failed", report);
    }
    return Status::OK();
}

// ============================================================================
// Compare two discovered_alleles maps for equivalence
// Returns true if they contain the same alleles with the same info
// ============================================================================
static bool discovered_alleles_equal(const discovered_alleles& a, const discovered_alleles& b) {
    if (a.size() != b.size()) return false;
    for (const auto& p : a) {
        auto it = b.find(p.first);
        if (it == b.end()) return false;
        if (p.second != it->second) return false;
    }
    return true;
}

// ============================================================================
// Incremental pipeline implementation
// ============================================================================

static std::shared_ptr<spdlog::logger> get_console() {
    auto logger = spdlog::get("GLnexus");
    if (!logger) {
        logger = spdlog::stderr_logger_mt("GLnexus_incremental");
    }
    return logger;
}
#define console get_console()

// Error-handling macro for the incremental CLI
#define HI(desc,expr) \
    s = expr; \
    if (s.bad()) { \
        console->error("incremental: failed to {}: {}", desc, s.str()); \
        return 1; \
    }

int incremental_steps(const vector<string>& new_gvcf_files,
                      const string& bedfilename,
                      const string& dbpath,
                      const string& config_name,
                      bool more_PL, bool squeeze, bool trim_uncalled_alleles,
                      const incremental_config& incr_cfg,
                      bool debug) {
    Status s;
    unifier_config unifier_cfg;
    genotyper_config genotyper_cfg;
    string cfg_txt, cfg_crc32c;

    size_t nr_threads = incr_cfg.nr_threads;
    size_t mem_budget = incr_cfg.mem_budget;

    if (nr_threads == 0) {
        nr_threads = std::thread::hardware_concurrency();
    }

    // Step 1: Load configuration
    HI("load unifier/genotyper configuration",
        cli::utils::load_config(console, config_name, unifier_cfg, genotyper_cfg, cfg_txt, cfg_crc32c,
                                more_PL, squeeze, trim_uncalled_alleles));

    // Step 2: Verify DB exists
    if (!cli::utils::check_dir_exists(dbpath)) {
        console->error("Database directory does not exist: {}. Use non-incremental mode for initial run.", dbpath);
        return 1;
    }

    // Phase 3: --compact-cache — clean up cached data and exit
    if (incr_cfg.compact_cache) {
        console->info("Compacting incremental cache...");
        unique_ptr<KeyValue::DB> db;
        RocksKeyValue::config opt;
        opt.mode = RocksKeyValue::OpenMode::NORMAL;
        opt.pfx = cli::utils::GLnexus_prefix_spec();
        opt.mem_budget = mem_budget;
        HI("open database", RocksKeyValue::Open(dbpath, opt, db));
        HI("ensure incremental collections", ensure_incremental_collections(db.get()));
        HI("clear genotype cache", clear_genotype_cache(db.get()));
        HI("flush database", db->flush());
        console->info("Cache compaction complete.");
        return 0;
    }

    // Step 3: Open DB briefly in NORMAL mode to read metadata, then close
    // (we'll need to close before db_bulk_load, which opens its own handle)
    vector<pair<string,size_t>> contigs;
    string old_sampleset;
    shared_ptr<const set<string>> old_samples;
    bool has_cache = false;
    {
        unique_ptr<KeyValue::DB> db;
        RocksKeyValue::config opt;
        opt.mode = RocksKeyValue::OpenMode::NORMAL;
        opt.pfx = cli::utils::GLnexus_prefix_spec();
        opt.mem_budget = mem_budget;
        opt.thread_budget = nr_threads;
        HI("open database for metadata", RocksKeyValue::Open(dbpath, opt, db));

        // Ensure incremental collections exist
        HI("ensure incremental collections", ensure_incremental_collections(db.get()));

        // Check config consistency
        {
            string cached_crc;
            Status cs = get_incremental_meta(db.get(), "config_crc32c", cached_crc);
            if (cs.ok() && cached_crc != cfg_crc32c) {
                if (incr_cfg.force_full) {
                    console->warn("Config CRC mismatch (cached={}, current={}). Forcing full re-run.",
                                  cached_crc, cfg_crc32c);
                } else {
                    console->error("Config CRC mismatch (cached={}, current={}). Use --force-full to override.",
                                   cached_crc, cfg_crc32c);
                    return 1;
                }
            }
        }

        // Get contigs
        {
            unique_ptr<BCFKeyValueData> data;
            HI("open BCFKeyValueData", BCFKeyValueData::Open(db.get(), data));
            HI("get contigs", data->contigs(contigs));
        }

        // Get old sampleset info (before importing new samples)
        {
            Status ms = get_incremental_meta(db.get(), "sampleset", old_sampleset);
            if (ms.ok() && !incr_cfg.force_full) {
                unique_ptr<BCFKeyValueData> data;
                HI("open BCFKeyValueData for old sampleset", BCFKeyValueData::Open(db.get(), data));
                Status ss = data->sampleset_samples(old_sampleset, old_samples);
                if (ss.ok()) {
                    has_cache = true;
                    console->info("Found cache from previous run ({} samples, sampleset={})",
                                  old_samples->size(), old_sampleset);
                } else {
                    console->warn("Could not load previous sampleset '{}': {}. Full pipeline.",
                                  old_sampleset, ss.str());
                }
            }
        }
        // DB handle released here
    }

    // Step 4: Import new gVCFs (db_bulk_load opens its own handle in BULK_LOAD mode)
    if (!new_gvcf_files.empty()) {
        vector<range> import_ranges; // empty = no range filter
        console->info("Importing {} new gVCF files...", new_gvcf_files.size());
        // Pass nullptr for db_out so db_bulk_load fully closes the handle (waits for compactions)
        HI("bulk load new gVCFs",
           cli::utils::db_bulk_load(console, mem_budget, nr_threads, new_gvcf_files, dbpath,
                                    import_ranges, contigs, nullptr, false));
    }

    // Step 5: Reopen DB in NORMAL mode for discovery/unification/caching
    unique_ptr<KeyValue::DB> db;
    {
        RocksKeyValue::config opt;
        opt.mode = RocksKeyValue::OpenMode::NORMAL;
        opt.pfx = cli::utils::GLnexus_prefix_spec();
        opt.mem_budget = mem_budget;
        opt.thread_budget = nr_threads;
        HI("reopen database", RocksKeyValue::Open(dbpath, opt, db));
    }
    // Re-ensure incremental collections (in case bulk_load recreated the DB)
    HI("ensure incremental collections", ensure_incremental_collections(db.get()));

    // Step 6: Get current sampleset and determine new samples
    string current_sampleset;
    shared_ptr<const set<string>> current_samples;
    unsigned sample_count = 0;
    set<string> added_samples;
    {
        unique_ptr<BCFKeyValueData> data;
        HI("open BCFKeyValueData for current sampleset", BCFKeyValueData::Open(db.get(), data));
        HI("get all_samples_sampleset", data->all_samples_sampleset(current_sampleset));
        HI("get current sampleset samples", data->sampleset_samples(current_sampleset, current_samples));
        sample_count = current_samples->size();
        console->info("Database now has {} samples (sampleset={})", sample_count, current_sampleset);

        if (has_cache && old_samples) {
            set_difference(current_samples->begin(), current_samples->end(),
                           old_samples->begin(), old_samples->end(),
                           inserter(added_samples, added_samples.end()));
            console->info("{} new samples detected since last incremental run", added_samples.size());
        }
    }

    // Step 10: Determine ranges
    vector<range> ranges;
    if (bedfilename.empty()) {
        console->warn("Processing full length of {} contigs, as no --bed was provided.", contigs.size());
        for (int rid = 0; rid < (int)contigs.size(); ++rid) {
            ranges.push_back(range(rid, 0, contigs[rid].second));
        }
    } else {
        HI("parse BED file", cli::utils::parse_bed_file(console, bedfilename, contigs, ranges));
    }

    // Step 11: Discovery
    // If we have cache AND new samples, do incremental discovery.
    // Otherwise, do full discovery and cache the results.
    vector<discovered_alleles> dsals_by_contig(contigs.size());
    vector<bool> contig_dirty(contigs.size(), false);

    if (has_cache && added_samples.empty()) {
        // ---- NO NEW SAMPLES: USE CACHED RESULTS ----
        console->info("No new samples detected. Using cached discovery results.");
        for (int rid = 0; rid < (int)contigs.size(); rid++) {
            Status cs = get_cached_discovered_alleles(db.get(), rid, dsals_by_contig[rid]);
            if (cs == StatusCode::NOT_FOUND) {
                // No cached alleles for this contig — that's normal (e.g., no variants)
            } else if (cs.bad()) {
                console->error("Failed to load cached discovered alleles for contig {}: {}",
                               contigs[rid].first, cs.str());
                return 1;
            }
            // All contigs are clean — unified sites cache is valid
        }
    } else if (has_cache && !added_samples.empty()) {
        // ---- INCREMENTAL DISCOVERY ----
        console->info("Running incremental discovery for {} new samples...", added_samples.size());

        // Create a temporary sampleset for the new samples
        string new_sampleset_name = "incr_new_" + to_string(sample_count);
        {
            unique_ptr<BCFKeyValueData> data;
            HI("open BCFKeyValueData for new sampleset", BCFKeyValueData::Open(db.get(), data));

            unique_ptr<MetadataCache> mcache;
            HI("create MetadataCache", MetadataCache::Start(*data, mcache));

            Status ns = data->new_sampleset(*mcache, new_sampleset_name, added_samples);
            if (ns == StatusCode::EXISTS) {
                console->info("Sampleset '{}' already exists, reusing.", new_sampleset_name);
            } else if (ns.bad()) {
                console->error("Failed to create new sampleset: {}", ns.str());
                return 1;
            }
        }

        // Discover alleles from the new-samples sampleset ONLY
        discovered_alleles new_dsals;
        unsigned new_N = 0;
        HI("discover alleles from new samples",
           cli::utils::discover_alleles_from_sampleset(
               console, nr_threads, db.get(), new_sampleset_name,
               ranges, contigs, new_dsals, new_N,
               unifier_cfg.min_allele_copy_number == 0));

        // Partition new discovered alleles by contig
        vector<discovered_alleles> new_dsals_by_contig(contigs.size());
        for (auto p = new_dsals.begin(); p != new_dsals.end(); new_dsals.erase(p++)) {
            UNPAIR(*p, al, dai);
            assert(al.pos.rid >= 0 && al.pos.rid < (int)contigs.size());
            new_dsals_by_contig[al.pos.rid][al] = dai;
        }

        // For each contig: load cached, merge, compare, cache updated
        for (int rid = 0; rid < (int)contigs.size(); rid++) {
            discovered_alleles cached;
            Status cs = get_cached_discovered_alleles(db.get(), rid, cached);

            if (cs == StatusCode::NOT_FOUND) {
                // No cache for this contig. Use new alleles if any.
                dsals_by_contig[rid] = move(new_dsals_by_contig[rid]);
                if (!dsals_by_contig[rid].empty()) {
                    contig_dirty[rid] = true;
                }
            } else if (cs.bad()) {
                console->error("Failed to load cached discovered alleles for contig {}: {}",
                               contigs[rid].first, cs.str());
                return 1;
            } else {
                // Merge cached + new
                discovered_alleles merged = cached;
                Status ms = merge_discovered_alleles(new_dsals_by_contig[rid], merged);
                if (ms.bad()) {
                    console->error("Failed to merge discovered alleles for contig {}: {}",
                                   contigs[rid].first, ms.str());
                    return 1;
                }

                // Compare old vs new
                if (!discovered_alleles_equal(cached, merged)) {
                    contig_dirty[rid] = true;
                    console->info("Contig {} has new/changed alleles (dirty)", contigs[rid].first);
                }
                dsals_by_contig[rid] = move(merged);
            }

            // Cache the updated discovered alleles
            if (!dsals_by_contig[rid].empty()) {
                HI("cache discovered alleles for " + contigs[rid].first,
                   put_cached_discovered_alleles(db.get(), rid, dsals_by_contig[rid]));
            }
        }

        size_t dirty_count = count(contig_dirty.begin(), contig_dirty.end(), true);
        console->info("Incremental discovery complete: {} of {} contigs are dirty",
                      dirty_count, contigs.size());
    } else {
        // ---- FULL DISCOVERY ----
        console->info("Running full discovery on all {} samples...", sample_count);

        discovered_alleles dsals;
        HI("discover alleles (full)",
           cli::utils::discover_alleles(console, nr_threads, db.get(), ranges, contigs, dsals,
                                        sample_count,
                                        unifier_cfg.min_allele_copy_number == 0));

        // Partition by contig
        for (auto p = dsals.begin(); p != dsals.end(); dsals.erase(p++)) {
            UNPAIR(*p, al, dai);
            assert(al.pos.rid >= 0 && al.pos.rid < (int)contigs.size());
            dsals_by_contig[al.pos.rid][al] = dai;
        }

        // Mark all contigs as dirty and cache
        for (int rid = 0; rid < (int)contigs.size(); rid++) {
            contig_dirty[rid] = true;
            if (!dsals_by_contig[rid].empty()) {
                HI("cache discovered alleles for " + contigs[rid].first,
                   put_cached_discovered_alleles(db.get(), rid, dsals_by_contig[rid]));
            }
        }
        console->info("Full discovery complete. All contigs marked dirty.");
    }

    // Step 12: Unification (parallel over contigs)
    console->info("Unifying sites...");
    auto nr_threads_m2 = nr_threads > 2 ? nr_threads - 2 : 1;
    ctpl::thread_pool unify_pool(nr_threads_m2);
    vector<future<Status>> unify_statuses;
    vector<vector<unified_site>> sites_by_contig(contigs.size());
    vector<unifier_stats> stats_by_contig(contigs.size());

    for (size_t i = 0; i < contigs.size(); i++) {
        unify_statuses.push_back(unify_pool.push([&, i](int tid) {
            if (!contig_dirty[i] && has_cache) {
                // Contig is clean — try to load cached unified sites
                Status cs = get_cached_unified_sites(db.get(), i, sites_by_contig[i]);
                if (cs.ok()) {
                    console->info("Contig {} is clean, using {} cached unified sites",
                                  contigs[i].first, sites_by_contig[i].size());
                    return Status::OK();
                }
                // If cache load fails, fall through to full unification
                console->warn("Failed to load cached unified sites for contig {}, re-unifying",
                              contigs[i].first);
            }

            // Full unification for this contig
            return cli::utils::unify_sites(console, unifier_cfg, contigs,
                                           dsals_by_contig[i], sample_count,
                                           sites_by_contig[i], stats_by_contig[i]);
        }));
    }

    // Collect unification results
    vector<unified_site> sites;
    unifier_stats stats;
    for (size_t i = 0; i < contigs.size(); i++) {
        HI("unify sites for " + contigs[i].first, unify_statuses[i].get());
        stats += stats_by_contig[i];

        // Cache unified sites for this contig (whether dirty or re-unified)
        if (contig_dirty[i] && !sites_by_contig[i].empty()) {
            HI("cache unified sites for " + contigs[i].first,
               put_cached_unified_sites(db.get(), i, sites_by_contig[i]));
        }

        auto& sites_i = sites_by_contig[i];
        sites.insert(sites.end(), make_move_iterator(sites_i.begin()),
                     make_move_iterator(sites_i.end()));
        sites_i.clear();
    }
    assert(std::is_sorted(sites.begin(), sites.end()));

    console->info("Unified to {} sites with {} ALT alleles. {} lost, {} filtered.",
                  sites.size(), stats.unified_alleles, stats.lost_alleles, stats.filtered_alleles);

    // Step 13: Update incremental metadata
    HI("save config CRC", put_incremental_meta(db.get(), "config_crc32c", cfg_crc32c));
    HI("save sampleset", put_incremental_meta(db.get(), "sampleset", current_sampleset));
    HI("save sample count", put_incremental_meta(db.get(), "sample_count", to_string(sample_count)));

    // Phase 2: Check if we can skip genotyping entirely (no new samples, all sites clean)
    bool any_dirty = false;
    for (size_t i = 0; i < contigs.size(); i++) {
        if (contig_dirty[i]) { any_dirty = true; break; }
    }

    // Step 14: Flush and close DB before genotyping (genotype() opens its own DB handle)
    console->info("Flushing database...");
    HI("flush database", db->flush());
    db.reset();

    if (has_cache && added_samples.empty() && !any_dirty && !incr_cfg.force_full) {
        // Phase 2: Nothing changed — output is identical to previous run
        console->info("No changes detected. Output would be identical to previous run. Skipping genotyping.");
        console->info("Incremental pipeline complete (no-op).");
        return 0;
    }

    // Step 15: Genotype all sites for all samples
    // Phase 3: Progress reporting
    console->info("Genotyping {} sites across {} samples...", sites.size(), sample_count);
    genotyper_cfg.output_residuals = debug;
    vector<string> hdr_lines = {
        ("##GLnexusConfigName=" + config_name),
        ("##GLnexusConfigCRC32C=" + cfg_crc32c),
        ("##GLnexusConfig=" + cfg_txt),
        "##GLnexusIncrementalMode=true"
    };
    string outfile("-");
    HI("genotype",
       cli::utils::genotype(console, mem_budget, nr_threads, dbpath, genotyper_cfg, sites, hdr_lines, outfile));

    // Phase 3: Validation mode — run full pipeline and compare
    if (incr_cfg.validate) {
        console->info("Validation mode: running full pipeline for comparison...");
        // Re-run full discovery on all samples
        {
            unique_ptr<KeyValue::DB> vdb;
            RocksKeyValue::config vopt;
            vopt.mode = RocksKeyValue::OpenMode::READ_ONLY;
            vopt.pfx = cli::utils::GLnexus_prefix_spec();
            vopt.mem_budget = mem_budget;
            vopt.thread_budget = nr_threads;
            HI("open database for validation", RocksKeyValue::Open(dbpath, vopt, vdb));

            discovered_alleles full_dsals;
            unsigned full_N = 0;
            HI("full discovery for validation",
               cli::utils::discover_alleles(console, nr_threads, vdb.get(), ranges, contigs,
                                            full_dsals, full_N, unifier_cfg.min_allele_copy_number == 0));

            // Compare allele counts
            size_t incr_total = 0;
            for (size_t i = 0; i < contigs.size(); i++) {
                incr_total += dsals_by_contig[i].size();
            }
            console->info("Validation: incremental discovered {} alleles, full discovered {} alleles",
                          incr_total, full_dsals.size());

            // Partition full alleles by contig and compare per-contig
            vector<discovered_alleles> full_by_contig(contigs.size());
            for (auto p = full_dsals.begin(); p != full_dsals.end(); full_dsals.erase(p++)) {
                UNPAIR(*p, al, dai);
                full_by_contig[al.pos.rid][al] = dai;
            }

            size_t mismatched_contigs = 0;
            for (size_t i = 0; i < contigs.size(); i++) {
                if (!discovered_alleles_equal(dsals_by_contig[i], full_by_contig[i])) {
                    console->warn("Validation: contig {} alleles differ (incremental={}, full={})",
                                  contigs[i].first, dsals_by_contig[i].size(), full_by_contig[i].size());
                    mismatched_contigs++;
                }
            }

            if (mismatched_contigs == 0) {
                console->info("Validation PASSED: incremental discovery matches full discovery");
            } else {
                console->error("Validation FAILED: {} contigs have mismatched alleles", mismatched_contigs);
            }
        }
    }

    console->info("Incremental pipeline complete.");
    return 0;
}

#undef HI
#undef console

}
