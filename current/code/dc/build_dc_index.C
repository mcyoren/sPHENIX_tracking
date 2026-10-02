#include "T_DigitalCurrent.C"

#include <TString.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

namespace dcindex
{
  struct FileRecord
  {
    int run = 0;
    int ebdc = -1;
    int stream = -1;
    int segment = -1;
    ULong64_t first_bco = 0;
    ULong64_t last_bco = 0;
    Long64_t entries = 0;
    std::string path;
  };

  bool valid_bco_entry(T_DigitalCurrent& t, Long64_t entry)
  {
    if (entry < 0 || entry >= t.GetEntries()) return false;
    if (t.LoadTree(entry) < 0) return false;
    t.GetEntry(entry);
    return t.dc_gtm_bco > 0 && t.dc_data_crc == t.dc_calc_crc;
  }

  bool find_first_last_bco(T_DigitalCurrent& t, ULong64_t& first, ULong64_t& last)
  {
    const Long64_t n = t.GetEntries();
    first = 0;
    last = 0;
    if (n <= 0) return false;

    for (Long64_t i = 0; i < n; ++i)
    {
      if (valid_bco_entry(t, i))
      {
        first = t.dc_gtm_bco;
        break;
      }
    }

    for (Long64_t i = n - 1; i >= 0; --i)
    {
      if (valid_bco_entry(t, i))
      {
        last = t.dc_gtm_bco;
        break;
      }
    }

    return first > 0 && last >= first;
  }
}

void build_dc_index(
    int run = 81558,
    const char* inputDir = "/sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run81558",
    const char* outputCsv = "dc_index_run81558.csv")
{
  using namespace dcindex;
  namespace fs = std::filesystem;

  const std::regex pattern(
      R"(^DCTree_DC-ebdc([0-9]{2})_([01])-([0-9]{8})-([0-9]{5})\.root$)");

  std::vector<FileRecord> records;

  if (!fs::exists(inputDir))
  {
    std::cerr << "ERROR: input directory does not exist: " << inputDir << std::endl;
    return;
  }

  for (const auto& entry : fs::directory_iterator(inputDir))
  {
    if (!entry.is_regular_file()) continue;

    const std::string name = entry.path().filename().string();
    std::smatch match;
    if (!std::regex_match(name, match, pattern)) continue;

    const int fileRun = std::stoi(match[3].str());
    if (fileRun != run) continue;

    FileRecord rec;
    rec.run = fileRun;
    rec.ebdc = std::stoi(match[1].str());
    rec.stream = std::stoi(match[2].str());
    rec.segment = std::stoi(match[4].str());
    rec.path = fs::absolute(entry.path()).string();

    std::cout << "Indexing ebdc=" << rec.ebdc
              << " stream=" << rec.stream
              << " segment=" << rec.segment
              << " : " << rec.path << std::endl;

    T_DigitalCurrent t(rec.path.c_str(), rec.ebdc);
    rec.entries = t.GetEntries();

    if (!find_first_last_bco(t, rec.first_bco, rec.last_bco))
    {
      std::cerr << "  WARNING: no valid nonzero CRC-good BCO found, skipping" << std::endl;
      continue;
    }

    std::cout << "  BCO [" << rec.first_bco << ", " << rec.last_bco
              << "] entries=" << rec.entries << std::endl;
    records.push_back(rec);
  }

  std::sort(records.begin(), records.end(), [](const FileRecord& a, const FileRecord& b)
  {
    if (a.ebdc != b.ebdc) return a.ebdc < b.ebdc;
    if (a.stream != b.stream) return a.stream < b.stream;
    if (a.first_bco != b.first_bco) return a.first_bco < b.first_bco;
    return a.segment < b.segment;
  });

  std::ofstream out(outputCsv);
  if (!out)
  {
    std::cerr << "ERROR: cannot write " << outputCsv << std::endl;
    return;
  }

  out << "run,ebdc,stream,segment,first_bco,last_bco,entries,path\n";
  for (const auto& rec : records)
  {
    out << rec.run << ','
        << rec.ebdc << ','
        << rec.stream << ','
        << rec.segment << ','
        << rec.first_bco << ','
        << rec.last_bco << ','
        << rec.entries << ','
        << rec.path << '\n';
  }
  out.close();

  bool present[24][2] = {{false}};
  for (const auto& rec : records)
  {
    if (rec.ebdc >= 0 && rec.ebdc < 24 && rec.stream >= 0 && rec.stream < 2)
      present[rec.ebdc][rec.stream] = true;
  }

  int missing = 0;
  for (int ebdc = 0; ebdc < 24; ++ebdc)
  {
    for (int stream = 0; stream < 2; ++stream)
    {
      if (!present[ebdc][stream])
      {
        std::cerr << "WARNING: no indexed files for ebdc " << ebdc
                  << " stream " << stream << std::endl;
        ++missing;
      }
    }
  }

  std::cout << "\nWrote " << records.size() << " indexed files to " << outputCsv << std::endl;
  if (missing == 0)
    std::cout << "All 48 EBDC/stream combinations are present." << std::endl;
  else
    std::cout << "Missing " << missing << " EBDC/stream combinations." << std::endl;
}
