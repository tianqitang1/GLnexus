#ifndef GLNEXUS_INCREMENTAL_H
#define GLNEXUS_INCREMENTAL_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include "types.h"
#include "data.h"
#include "KeyValue.h"
#include "BCFKeyValueData.h"
#include "unifier.h"
#include "spdlog/spdlog.h"

namespace GLnexus {

// ============================================================================
// Binary serialization for discovered_alleles
// ============================================================================

Status serialize_discovered_alleles(const discovered_alleles& dsals, std::string& out);
Status deserialize_discovered_alleles(const char* data, size_t size, discovered_alleles& ans);
Status deserialize_discovered_alleles(const std::string& data, discovered_alleles& ans);

// ============================================================================
// Binary serialization for unified_site
// ============================================================================

Status serialize_unified_sites(const std::vector<unified_site>& sites, std::string& out);
Status deserialize_unified_sites(const char* data, size_t size, std::vector<unified_site>& ans);
Status deserialize_unified_sites(const std::string& data, std::vector<unified_site>& ans);

// ============================================================================
// Incremental column family helpers
// ============================================================================

extern const char* DISCOVERED_CF;
extern const char* UNIFIED_CF;
extern const char* GENOTYPE_CACHE_CF;   // Phase 2: per-site BCF record cache
extern const char* INCREMENTAL_META_CF;

Status ensure_incremental_collections(KeyValue::DB* db);
std::string contig_key(int rid);

// ============================================================================
// Incremental pipeline metadata
// ============================================================================

Status put_incremental_meta(KeyValue::DB* db, const std::string& key, const std::string& value);
Status get_incremental_meta(KeyValue::DB* db, const std::string& key, std::string& value);

// ============================================================================
// Cached discovered_alleles (per contig)
// ============================================================================

Status put_cached_discovered_alleles(KeyValue::DB* db, int rid, const discovered_alleles& dsals);
Status get_cached_discovered_alleles(KeyValue::DB* db, int rid, discovered_alleles& dsals);

// ============================================================================
// Cached unified_sites (per contig)
// ============================================================================

Status put_cached_unified_sites(KeyValue::DB* db, int rid, const std::vector<unified_site>& sites);
Status get_cached_unified_sites(KeyValue::DB* db, int rid, std::vector<unified_site>& sites);

// ============================================================================
// Phase 2: Genotype cache (per site)
// ============================================================================

// Encode a site key from rid + position for the genotype cache
std::string site_key(int rid, int beg, int end);

// Store a genotyped BCF record for a site
Status put_cached_genotype(KeyValue::DB* db, const range& pos, const std::string& bcf_data);
Status get_cached_genotype(KeyValue::DB* db, const range& pos, std::string& bcf_data);

// Invalidate all cached genotypes (e.g., after config change)
Status clear_genotype_cache(KeyValue::DB* db);

// ============================================================================
// Phase 3: Validation - compare two BCF files at the genotype level
// ============================================================================

// Compare genotypes from two BCF files field by field (GT, GQ, DP).
// Returns OK if they match, Invalid with description if they differ.
Status validate_bcf_genotypes(const std::string& file_a, const std::string& file_b,
                              std::string& report);

// ============================================================================
// Incremental pipeline configuration
// ============================================================================

struct incremental_config {
    size_t mem_budget = 0;
    size_t nr_threads = 0;
    bool force_full = false;     // force full re-run ignoring cache
    bool validate = false;       // Phase 3: compare incremental vs full output
    bool compact_cache = false;  // Phase 3: compact/clean up cache data
};

// Run the incremental pipeline.
// Returns 0 on success, 1 on failure.
int incremental_steps(const std::vector<std::string>& new_gvcf_files,
                      const std::string& bedfilename,
                      const std::string& dbpath,
                      const std::string& config_name,
                      bool more_PL, bool squeeze, bool trim_uncalled_alleles,
                      const incremental_config& incr_cfg,
                      bool debug);

}

#endif
