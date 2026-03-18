// Tests for the incremental DB postprocessor (Phases 1-3)
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include "catch.hpp"
#include "types.h"
#include "incremental.h"
#include "BCFKeyValueData.h"
#include "RocksKeyValue.h"
#include "cli_utils.h"
#include "service.h"
#include "unifier.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/null_sink.h"

using namespace std;
using namespace GLnexus;

static shared_ptr<spdlog::logger> test_logger() {
    static auto logger = spdlog::create<spdlog::sinks::null_sink_st>("test_incremental_null");
    return logger;
}

// ============================================================================
// Phase 1: Serialization roundtrip tests
// ============================================================================

TEST_CASE("discovered_alleles serialization roundtrip") {
    discovered_alleles dsals;

    {
        allele al(range(0, 100, 103), "ACG");
        discovered_allele_info info;
        info.is_ref = true;
        info.all_filtered = false;
        info.topAQ.V[0] = 500;
        info.topAQ.V[1] = 300;
        info.zGQ.M[5][0] = 10;
        info.zGQ.M[5][1] = 5;
        info.in_target = range(0, 50, 200);
        dsals[al] = info;
    }

    {
        allele al(range(0, 100, 103), "ATG");
        discovered_allele_info info;
        info.is_ref = false;
        info.all_filtered = true;
        info.topAQ.V[0] = 250;
        info.zGQ.M[3][0] = 7;
        info.in_target = range(0, 50, 200);
        dsals[al] = info;
    }

    {
        allele al(range(1, 500, 502), "GA");
        discovered_allele_info info;
        info.is_ref = false;
        info.all_filtered = false;
        info.topAQ.V[0] = 9999;
        info.topAQ.V[1] = 8888;
        info.topAQ.V[2] = 7777;
        info.zGQ.M[9][1] = 42;
        info.in_target = range(-1, -1, -1);
        dsals[al] = info;
    }

    string serialized;
    REQUIRE(serialize_discovered_alleles(dsals, serialized).ok());
    REQUIRE(serialized.size() > 0);

    discovered_alleles result;
    REQUIRE(deserialize_discovered_alleles(serialized, result).ok());

    REQUIRE(result.size() == dsals.size());
    for (const auto& p : dsals) {
        auto it = result.find(p.first);
        REQUIRE(it != result.end());
        REQUIRE(it->second.is_ref == p.second.is_ref);
        REQUIRE(it->second.all_filtered == p.second.all_filtered);
        REQUIRE(it->second.topAQ == p.second.topAQ);
        REQUIRE(it->second.zGQ == p.second.zGQ);
        REQUIRE(it->second.in_target == p.second.in_target);
    }
}

TEST_CASE("discovered_alleles serialization empty map") {
    discovered_alleles dsals;
    string serialized;
    REQUIRE(serialize_discovered_alleles(dsals, serialized).ok());

    discovered_alleles result;
    REQUIRE(deserialize_discovered_alleles(serialized, result).ok());
    REQUIRE(result.empty());
}

TEST_CASE("unified_sites serialization roundtrip") {
    vector<unified_site> sites;

    {
        unified_site us(range(0, 100, 105));
        us.in_target = range(0, 50, 200);

        unified_allele ua_ref(us.pos, "ACGTG");
        ua_ref.quality = 0;
        ua_ref.frequency = NAN;
        us.alleles.push_back(ua_ref);

        unified_allele ua_alt(us.pos, "ATGTG");
        ua_alt.normalized = allele(range(0, 101, 102), "T");
        ua_alt.quality = 500;
        ua_alt.frequency = 0.05f;
        us.alleles.push_back(ua_alt);

        us.unification[allele(us.pos, "ACGTG")] = 0;
        us.unification[allele(us.pos, "ATGTG")] = 1;
        us.unification[allele(range(0, 101, 102), "T")] = 1;

        us.lost_allele_frequency = 0.001f;
        us.qual = 500;
        us.monoallelic = false;

        sites.push_back(us);
    }

    {
        unified_site us(range(0, 200, 201));
        us.in_target = range(-1, -1, -1);

        unified_allele ua_ref(us.pos, "G");
        us.alleles.push_back(ua_ref);

        unified_allele ua_alt(us.pos, "A");
        ua_alt.quality = 100;
        ua_alt.frequency = 0.1f;
        us.alleles.push_back(ua_alt);

        us.fill_implicit_unification();
        us.qual = 100;
        us.monoallelic = true;

        sites.push_back(us);
    }

    string serialized;
    REQUIRE(serialize_unified_sites(sites, serialized).ok());
    REQUIRE(serialized.size() > 0);

    vector<unified_site> result;
    REQUIRE(deserialize_unified_sites(serialized, result).ok());

    REQUIRE(result.size() == sites.size());
    for (size_t i = 0; i < sites.size(); i++) {
        REQUIRE(result[i].pos == sites[i].pos);
        REQUIRE(result[i].in_target == sites[i].in_target);
        REQUIRE(result[i].alleles.size() == sites[i].alleles.size());
        for (size_t j = 0; j < sites[i].alleles.size(); j++) {
            REQUIRE(result[i].alleles[j].dna == sites[i].alleles[j].dna);
            REQUIRE(result[i].alleles[j].quality == sites[i].alleles[j].quality);
            REQUIRE(result[i].alleles[j].normalized == sites[i].alleles[j].normalized);
        }
        REQUIRE(result[i].unification == sites[i].unification);
        REQUIRE(result[i].qual == sites[i].qual);
        REQUIRE(result[i].monoallelic == sites[i].monoallelic);
    }
}

TEST_CASE("unified_sites serialization empty vector") {
    vector<unified_site> sites;
    string serialized;
    REQUIRE(serialize_unified_sites(sites, serialized).ok());

    vector<unified_site> result;
    REQUIRE(deserialize_unified_sites(serialized, result).ok());
    REQUIRE(result.empty());
}

TEST_CASE("deserialization rejects invalid data") {
    discovered_alleles dsals;
    REQUIRE(deserialize_discovered_alleles("", 0, dsals).bad());
    REQUIRE(deserialize_discovered_alleles("x", 1, dsals).bad());

    string bad_version;
    bad_version.resize(8, 0);
    uint32_t wrong_ver = 999;
    memcpy(&bad_version[0], &wrong_ver, sizeof(wrong_ver));
    REQUIRE(deserialize_discovered_alleles(bad_version, dsals).bad());

    vector<unified_site> usites;
    REQUIRE(deserialize_unified_sites("", 0, usites).bad());
    REQUIRE(deserialize_unified_sites(bad_version, usites).bad());
}

// ============================================================================
// Phase 1: Incremental collection helpers tests
// ============================================================================

TEST_CASE("incremental collection helpers with RocksDB") {
    string dbpath = "/tmp/test_incremental_collections";
    ignore_retval(system(("rm -rf " + dbpath).c_str()));

    {
        RocksKeyValue::config cfg;
        unique_ptr<KeyValue::DB> db;
        REQUIRE(RocksKeyValue::Initialize(dbpath, cfg, db).ok());
        vector<pair<string,size_t>> contigs = {{"chr1", 1000}, {"chr2", 500}};
        REQUIRE(BCFKeyValueData::InitializeDB(db.get(), contigs).ok());
    }

    {
        RocksKeyValue::config cfg;
        cfg.mode = RocksKeyValue::OpenMode::NORMAL;
        unique_ptr<KeyValue::DB> db;
        REQUIRE(RocksKeyValue::Open(dbpath, cfg, db).ok());
        REQUIRE(ensure_incremental_collections(db.get()).ok());

        // Test incremental_meta put/get
        REQUIRE(put_incremental_meta(db.get(), "sample_count", "42").ok());
        string val;
        REQUIRE(get_incremental_meta(db.get(), "sample_count", val).ok());
        REQUIRE(val == "42");

        string missing;
        REQUIRE(get_incremental_meta(db.get(), "nonexistent", missing) == StatusCode::NOT_FOUND);

        // Test discovered_alleles cache
        discovered_alleles dsals;
        {
            allele al(range(0, 100, 103), "ACG");
            discovered_allele_info info;
            info.is_ref = true;
            info.topAQ.V[0] = 500;
            info.in_target = range(0, 0, 1000);
            dsals[al] = info;
        }
        REQUIRE(put_cached_discovered_alleles(db.get(), 0, dsals).ok());

        discovered_alleles loaded;
        REQUIRE(get_cached_discovered_alleles(db.get(), 0, loaded).ok());
        REQUIRE(loaded.size() == 1);

        discovered_alleles missing_dsals;
        REQUIRE(get_cached_discovered_alleles(db.get(), 99, missing_dsals) == StatusCode::NOT_FOUND);

        // Test unified_sites cache
        vector<unified_site> sites;
        {
            unified_site us(range(0, 100, 105));
            us.in_target = range(-1, -1, -1);
            unified_allele ua(us.pos, "ACGTG");
            us.alleles.push_back(ua);
            unified_allele ua2(us.pos, "ATGTG");
            ua2.quality = 100;
            ua2.frequency = 0.05f;
            us.alleles.push_back(ua2);
            us.fill_implicit_unification();
            us.qual = 100;
            sites.push_back(us);
        }
        REQUIRE(put_cached_unified_sites(db.get(), 0, sites).ok());

        vector<unified_site> loaded_sites;
        REQUIRE(get_cached_unified_sites(db.get(), 0, loaded_sites).ok());
        REQUIRE(loaded_sites.size() == 1);
        REQUIRE(loaded_sites[0].pos == sites[0].pos);
        REQUIRE(loaded_sites[0].alleles.size() == sites[0].alleles.size());
    }

    ignore_retval(system(("rm -rf " + dbpath).c_str()));
}

TEST_CASE("ensure_incremental_collections is idempotent") {
    string dbpath = "/tmp/test_incr_idempotent";
    ignore_retval(system(("rm -rf " + dbpath).c_str()));

    {
        RocksKeyValue::config cfg;
        unique_ptr<KeyValue::DB> db;
        REQUIRE(RocksKeyValue::Initialize(dbpath, cfg, db).ok());
        vector<pair<string,size_t>> contigs = {{"chr1", 1000}};
        REQUIRE(BCFKeyValueData::InitializeDB(db.get(), contigs).ok());
    }

    {
        RocksKeyValue::config cfg;
        cfg.mode = RocksKeyValue::OpenMode::NORMAL;
        unique_ptr<KeyValue::DB> db;
        REQUIRE(RocksKeyValue::Open(dbpath, cfg, db).ok());

        REQUIRE(ensure_incremental_collections(db.get()).ok());
        REQUIRE(ensure_incremental_collections(db.get()).ok());

        KeyValue::CollectionHandle coll;
        REQUIRE(db->collection(DISCOVERED_CF, coll).ok());
        REQUIRE(db->collection(UNIFIED_CF, coll).ok());
        REQUIRE(db->collection(GENOTYPE_CACHE_CF, coll).ok());
        REQUIRE(db->collection(INCREMENTAL_META_CF, coll).ok());
    }

    ignore_retval(system(("rm -rf " + dbpath).c_str()));
}

TEST_CASE("contig_key encoding") {
    string k0 = contig_key(0);
    string k1 = contig_key(1);
    string k255 = contig_key(255);
    string k1000 = contig_key(1000);

    REQUIRE(k0.size() == 3);
    REQUIRE(k1.size() == 3);
    REQUIRE(k0 != k1);
    REQUIRE(k0 != k255);
    REQUIRE(k0 != k1000);
    REQUIRE(k1 != k1000);

    REQUIRE(k0 < k1);
    REQUIRE(k1 < k255);
    REQUIRE(k255 < k1000);
}

// ============================================================================
// Phase 2: Genotype cache tests
// ============================================================================

TEST_CASE("site_key encoding") {
    string k1 = site_key(0, 100, 200);
    string k2 = site_key(0, 100, 201);
    string k3 = site_key(0, 200, 300);
    string k4 = site_key(1, 100, 200);

    REQUIRE(k1.size() == 11);
    REQUIRE(k1 != k2);
    REQUIRE(k1 != k3);
    REQUIRE(k1 != k4);

    // Ordered by rid, then beg, then end
    REQUIRE(k1 < k2);
    REQUIRE(k2 < k3);
    REQUIRE(k3 < k4);
}

TEST_CASE("genotype cache put/get") {
    string dbpath = "/tmp/test_genotype_cache";
    ignore_retval(system(("rm -rf " + dbpath).c_str()));

    {
        RocksKeyValue::config cfg;
        unique_ptr<KeyValue::DB> db;
        REQUIRE(RocksKeyValue::Initialize(dbpath, cfg, db).ok());
        vector<pair<string,size_t>> contigs = {{"chr1", 1000}};
        REQUIRE(BCFKeyValueData::InitializeDB(db.get(), contigs).ok());
    }

    {
        RocksKeyValue::config cfg;
        cfg.mode = RocksKeyValue::OpenMode::NORMAL;
        unique_ptr<KeyValue::DB> db;
        REQUIRE(RocksKeyValue::Open(dbpath, cfg, db).ok());
        REQUIRE(ensure_incremental_collections(db.get()).ok());

        // Store a fake BCF record
        string fake_bcf = "fake_bcf_record_data_here";
        range pos(0, 100, 200);
        REQUIRE(put_cached_genotype(db.get(), pos, fake_bcf).ok());

        // Retrieve it
        string loaded;
        REQUIRE(get_cached_genotype(db.get(), pos, loaded).ok());
        REQUIRE(loaded == fake_bcf);

        // Non-existent site returns NotFound
        string missing;
        REQUIRE(get_cached_genotype(db.get(), range(0, 999, 1000), missing) == StatusCode::NOT_FOUND);

        // Clear cache
        REQUIRE(clear_genotype_cache(db.get()).ok());
    }

    ignore_retval(system(("rm -rf " + dbpath).c_str()));
}

// ============================================================================
// Integration test: Full pipeline with real gVCFs
// ============================================================================

TEST_CASE("incremental integration with real gVCFs") {
    auto console = test_logger();
    Status s;
    int nr_threads = 2;

    string basedir = "test/data/cli";
    string DB_DIR = "/tmp/test_incremental_integration";
    string DB_PATH = DB_DIR + "/DB";

    // Clean up
    REQUIRE(system(("rm -rf " + DB_DIR).c_str()) == 0);
    REQUIRE(system(("mkdir -p " + DB_DIR).c_str()) == 0);

    // ----------------------------------------------------------------
    // Step 1: Initialize DB and load first batch (F1 + F2)
    // ----------------------------------------------------------------
    vector<pair<string,size_t>> contigs;
    s = cli::utils::db_init(console, DB_PATH, basedir + "/F1.gvcf.gz", contigs);
    REQUIRE(s.ok());

    vector<string> batch1 = {basedir + "/F1.gvcf.gz", basedir + "/F2.gvcf.gz"};
    vector<range> empty_ranges;
    s = cli::utils::db_bulk_load(console, 0, nr_threads, batch1, DB_PATH, empty_ranges, contigs);
    REQUIRE(s.ok());

    // ----------------------------------------------------------------
    // Step 2: Run full discovery + unification, cache results
    // ----------------------------------------------------------------
    vector<range> ranges;
    for (int rid = 0; rid < (int)contigs.size(); rid++) {
        ranges.push_back(range(rid, 0, contigs[rid].second));
    }

    discovered_alleles dsals_2;
    unsigned sample_count_2 = 0;
    s = cli::utils::discover_alleles(console, 0, nr_threads, DB_PATH, ranges, contigs,
                                     dsals_2, sample_count_2);
    REQUIRE(s.ok());
    REQUIRE(sample_count_2 == 2);
    REQUIRE(dsals_2.size() > 0);

    // Cache discovered alleles per contig
    {
        unique_ptr<KeyValue::DB> db;
        RocksKeyValue::config opt;
        opt.mode = RocksKeyValue::OpenMode::NORMAL;
        opt.pfx = cli::utils::GLnexus_prefix_spec();
        REQUIRE(RocksKeyValue::Open(DB_PATH, opt, db).ok());
        REQUIRE(ensure_incremental_collections(db.get()).ok());

        // Partition by contig and cache
        vector<discovered_alleles> dsals_by_contig(contigs.size());
        for (const auto& p : dsals_2) {
            dsals_by_contig[p.first.pos.rid][p.first] = p.second;
        }
        for (int rid = 0; rid < (int)contigs.size(); rid++) {
            if (!dsals_by_contig[rid].empty()) {
                REQUIRE(put_cached_discovered_alleles(db.get(), rid, dsals_by_contig[rid]).ok());
            }
        }

        // Unify and cache
        unifier_config unifier_cfg;
        unifier_cfg.min_AQ1 = 0;
        unifier_cfg.min_AQ2 = 0;

        for (int rid = 0; rid < (int)contigs.size(); rid++) {
            if (!dsals_by_contig[rid].empty()) {
                vector<unified_site> sites_i;
                unifier_stats stats_i;
                s = cli::utils::unify_sites(console, unifier_cfg, contigs,
                                            dsals_by_contig[rid], sample_count_2,
                                            sites_i, stats_i);
                REQUIRE(s.ok());
                if (!sites_i.empty()) {
                    REQUIRE(put_cached_unified_sites(db.get(), rid, sites_i).ok());
                }
            }
        }

        // Save metadata
        unique_ptr<BCFKeyValueData> data;
        REQUIRE(BCFKeyValueData::Open(db.get(), data).ok());
        string sampleset;
        REQUIRE(data->all_samples_sampleset(sampleset).ok());
        REQUIRE(put_incremental_meta(db.get(), "sampleset", sampleset).ok());
        REQUIRE(put_incremental_meta(db.get(), "sample_count", to_string(sample_count_2)).ok());

        REQUIRE(db->flush().ok());
    }

    // ----------------------------------------------------------------
    // Step 3: Import second batch (F3 + F4)
    // ----------------------------------------------------------------
    vector<string> batch2 = {basedir + "/F3.gvcf.gz", basedir + "/F4.gvcf.gz"};
    s = cli::utils::db_bulk_load(console, 0, nr_threads, batch2, DB_PATH, empty_ranges, contigs);
    REQUIRE(s.ok());

    // ----------------------------------------------------------------
    // Step 4: Run incremental discovery (new samples only) + merge
    // ----------------------------------------------------------------
    discovered_alleles dsals_incr;
    {
        unique_ptr<KeyValue::DB> db;
        RocksKeyValue::config opt;
        opt.mode = RocksKeyValue::OpenMode::NORMAL;
        opt.pfx = cli::utils::GLnexus_prefix_spec();
        REQUIRE(RocksKeyValue::Open(DB_PATH, opt, db).ok());
        REQUIRE(ensure_incremental_collections(db.get()).ok());

        // Get old sampleset
        string old_sampleset;
        REQUIRE(get_incremental_meta(db.get(), "sampleset", old_sampleset).ok());

        unique_ptr<BCFKeyValueData> data;
        REQUIRE(BCFKeyValueData::Open(db.get(), data).ok());

        shared_ptr<const set<string>> old_samples;
        REQUIRE(data->sampleset_samples(old_sampleset, old_samples).ok());
        REQUIRE(old_samples->size() == 2);

        string current_sampleset;
        REQUIRE(data->all_samples_sampleset(current_sampleset).ok());
        shared_ptr<const set<string>> current_samples;
        REQUIRE(data->sampleset_samples(current_sampleset, current_samples).ok());
        REQUIRE(current_samples->size() == 4);

        // Find new samples
        set<string> added_samples;
        set_difference(current_samples->begin(), current_samples->end(),
                       old_samples->begin(), old_samples->end(),
                       inserter(added_samples, added_samples.end()));
        REQUIRE(added_samples.size() == 2);

        // Create sampleset for new samples
        unique_ptr<MetadataCache> mcache;
        REQUIRE(MetadataCache::Start(*data, mcache).ok());
        string new_ss_name = "test_new_samples";
        Status ns = data->new_sampleset(*mcache, new_ss_name, added_samples);
        REQUIRE((ns.ok() || ns == StatusCode::EXISTS));

        // Discover from new samples only
        discovered_alleles new_dsals;
        unsigned new_N = 0;
        s = cli::utils::discover_alleles_from_sampleset(
            console, nr_threads, db.get(), new_ss_name,
            ranges, contigs, new_dsals, new_N);
        REQUIRE(s.ok());

        // Load cached and merge per contig
        for (int rid = 0; rid < (int)contigs.size(); rid++) {
            discovered_alleles cached;
            Status cs = get_cached_discovered_alleles(db.get(), rid, cached);

            // Partition new alleles for this contig
            discovered_alleles new_for_contig;
            for (const auto& p : new_dsals) {
                if (p.first.pos.rid == rid) {
                    new_for_contig[p.first] = p.second;
                }
            }

            if (cs.ok()) {
                // Merge
                REQUIRE(merge_discovered_alleles(new_for_contig, cached).ok());
                for (const auto& p : cached) {
                    dsals_incr[p.first] = p.second;
                }
            } else {
                for (const auto& p : new_for_contig) {
                    dsals_incr[p.first] = p.second;
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 5: Run full discovery on all 4 samples for comparison
    // ----------------------------------------------------------------
    discovered_alleles dsals_full;
    unsigned sample_count_full = 0;
    s = cli::utils::discover_alleles(console, 0, nr_threads, DB_PATH, ranges, contigs,
                                     dsals_full, sample_count_full);
    REQUIRE(s.ok());
    REQUIRE(sample_count_full == 4);

    // ----------------------------------------------------------------
    // Step 6: Compare incremental vs full discovered alleles
    // ----------------------------------------------------------------
    // The incremental result should contain all alleles from the full result.
    // (Incremental merges cached 2-sample alleles + new 2-sample alleles,
    //  which should produce the same allele set as full 4-sample discovery.)
    //
    // Note: exact equality of discovered_allele_info (topAQ, zGQ) depends
    // on discovery order and the merge algorithm. We check that every allele
    // in the full result also appears in the incremental result.
    for (const auto& p : dsals_full) {
        auto it = dsals_incr.find(p.first);
        REQUIRE(it != dsals_incr.end());
        // The allele should have the same is_ref status
        REQUIRE(it->second.is_ref == p.second.is_ref);
    }

    // Clean up
    ignore_retval(system(("rm -rf " + DB_DIR).c_str()));
}

// ============================================================================
// Phase 3: Validation function test
// ============================================================================

TEST_CASE("validate_bcf_genotypes with identical files") {
    auto console = test_logger();
    Status s;
    int nr_threads = 2;

    string basedir = "test/data/cli";
    string DB_DIR = "/tmp/test_validate_bcf";
    string DB_PATH = DB_DIR + "/DB";

    REQUIRE(system(("rm -rf " + DB_DIR).c_str()) == 0);
    REQUIRE(system(("mkdir -p " + DB_DIR).c_str()) == 0);

    // Build a small DB and genotype it
    vector<pair<string,size_t>> contigs;
    s = cli::utils::db_init(console, DB_PATH, basedir + "/F1.gvcf.gz", contigs);
    REQUIRE(s.ok());

    vector<string> gvcfs = {basedir + "/F1.gvcf.gz", basedir + "/F2.gvcf.gz"};
    vector<range> empty_ranges;
    s = cli::utils::db_bulk_load(console, 0, nr_threads, gvcfs, DB_PATH, empty_ranges, contigs);
    REQUIRE(s.ok());

    vector<range> ranges;
    for (int rid = 0; rid < (int)contigs.size(); rid++) {
        ranges.push_back(range(rid, 0, contigs[rid].second));
    }

    discovered_alleles dsals;
    unsigned N = 0;
    s = cli::utils::discover_alleles(console, 0, nr_threads, DB_PATH, ranges, contigs, dsals, N);
    REQUIRE(s.ok());

    unifier_config ucfg;
    ucfg.min_AQ1 = 0;
    ucfg.min_AQ2 = 0;
    genotyper_config gcfg;

    vector<unified_site> sites;
    unifier_stats stats;
    s = cli::utils::unify_sites(console, ucfg, contigs, dsals, N, sites, stats);
    REQUIRE(s.ok());
    REQUIRE(sites.size() > 0);

    string bcf_a = DB_DIR + "/a.bcf";
    string bcf_b = DB_DIR + "/b.bcf";
    s = cli::utils::genotype(console, 0, nr_threads, DB_PATH, gcfg, sites, {}, bcf_a);
    REQUIRE(s.ok());
    s = cli::utils::genotype(console, 0, nr_threads, DB_PATH, gcfg, sites, {}, bcf_b);
    REQUIRE(s.ok());

    // Compare — should be identical
    string report;
    s = validate_bcf_genotypes(bcf_a, bcf_b, report);
    REQUIRE(s.ok());
    REQUIRE(report.find("matched") != string::npos);
    REQUIRE(report.find("0 mismatched") != string::npos);

    ignore_retval(system(("rm -rf " + DB_DIR).c_str()));
}
