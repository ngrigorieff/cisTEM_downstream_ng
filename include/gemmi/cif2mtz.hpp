//The contents of this file are covered by the Mozilla Public License v2, a copy of which is included in include/LICENSE_MOZILLAv2.txt
// Copyright 2021 Global Phasing Ltd.
//
// A class for converting SF-mmCIF to MTZ (merged or unmerged).

#ifndef GEMMI_CIF2MTZ_HPP_
#define GEMMI_CIF2MTZ_HPP_

#include <ostream>
#include <map>
#include <set>
#include <utility>
#include "cifdoc.hpp"   // for Loop, as_int, ...
#include "fail.hpp"     // for fail
#include "mtz.hpp"      // for Mtz
#include "numb.hpp"     // for as_number
#include "refln.hpp"    // for ReflnBlock
#include "version.hpp"  // for GEMMI_VERSION

namespace gemmi {

struct CifToMtz {
  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "pdbx_r_free_flag FreeR_flag I 0",
      "status FreeR_flag s 0", // s is a special flag
      "intensity_meas IMEAN J 1",
      "intensity_sigma SIGIMEAN Q 1",
      "pdbx_I_plus I(+) K 1",
      "pdbx_I_plus_sigma SIGI(+) M 1",
      "pdbx_I_minus I(-) K 1",
      "pdbx_I_minus_sigma SIGI(-) M 1",
      "F_meas FP F 1",
      "F_meas_au FP F 1",
      "F_meas_sigma SIGFP Q 1",
      "F_meas_sigma_au SIGFP Q 1",
      "pdbx_F_plus F(+) G 1",
      "pdbx_F_plus_sigma SIGF(+) L 1",
      "pdbx_F_minus F(-) G 1",
      "pdbx_F_minus_sigma SIGF(-) L 1",
      "pdbx_anom_difference DP D 1",
      "pdbx_anom_difference_sigma SIGDP Q 1",
      "F_calc FC F 1",
      "F_calc_au FC F 1",
      "phase_calc PHIC P 1",
      "fom FOM W 1",
      "weight FOM W 1",
      "pdbx_HL_A_iso HLA A 1",
      "pdbx_HL_B_iso HLB A 1",
      "pdbx_HL_C_iso HLC A 1",
      "pdbx_HL_D_iso HLD A 1",
      "pdbx_FWT FWT F 1",
      "pdbx_PHWT PHWT P 1",
      "pdbx_DELFWT DELFWT F 1",
      "pdbx_DELPHWT DELPHWT P 1",
      nullptr
    };
    static const char* unmerged[] = {
      "intensity_meas I J 1",  // for unmerged data is category refln
      "intensity_net I J 1",
      "intensity_sigma SIGI Q 1",
      "pdbx_detector_x XDET R 1",
      "pdbx_detector_y YDET R 1",
      "pdbx_scan_angle ROT R 1",
      nullptr
    };
    return for_merged ? merged : unmerged;
  };

  struct Entry {
    std::string refln_tag;
    std::string col_label;
    char col_type;
    int dataset_id;

    Entry(const std::string& line) {
      std::vector<std::string> tokens;
      tokens.reserve(4);
      split_str_into_multi(line, " \t\r\n", tokens);
      if (tokens.size() != 4)
        fail("line should have 4 words: " + line);
      if (tokens[2].size() != 1 || tokens[3].size() != 1 ||
          (tokens[3][0] != '0' && tokens[3][0] != '1'))
        fail("incorrect line: " + line);
      refln_tag = tokens[0];
      col_label = tokens[1];
      col_type = tokens[2][0];
      dataset_id = tokens[3][0] - '0';
    }
  };

  // Alternative mmCIF tags for the same MTZ label should be consecutive
  bool verbose = false;
  bool force_unmerged = false;
  std::string title;
  std::vector<std::string> history = { "From gemmi-cif2mtz " GEMMI_VERSION };
  std::vector<std::string> spec_lines;

  Mtz convert_block_to_mtz(const ReflnBlock& rb, std::ostream& out) const {
    Mtz mtz;
    mtz.title = title.empty() ? "Converted from mmCIF block " + rb.block.name : title;
    if (!history.empty()) {
      mtz.history.reserve(mtz.history.size() + history.size());
      mtz.history.insert(mtz.history.end(), history.begin(), history.end());
    }
    mtz.cell = rb.cell;
    mtz.spacegroup = rb.spacegroup;
    mtz.add_dataset("HKL_base");
    mtz.add_dataset("unknown").wavelength = rb.wavelength;
    const cif::Loop* loop = rb.refln_loop ? rb.refln_loop : rb.diffrn_refln_loop;
    if (!loop)
      fail("_refln category not found in mmCIF block: " + rb.block.name);
    if (verbose)
      out << "Searching tags with known MTZ equivalents ...\n";
    bool uses_status = false;
    std::vector<int> indices;
    std::string tag = loop->tags[0];
    const size_t len = tag.find('.') + 1;

    bool unmerged = force_unmerged || !rb.refln_loop;

    std::vector<Entry> spec_entries;
    if (!spec_lines.empty()) {
      spec_entries.reserve(spec_lines.size());
      for (const std::string& line : spec_lines)
        spec_entries.emplace_back(line);
    } else {
      const char** line = default_spec(!unmerged);
      for (; *line != nullptr; ++line)
        spec_entries.emplace_back(*line);
    }

    // always start with H, K, L
    tag.replace(len, std::string::npos, "index_h");
    for (char c : {'h', 'k', 'l'}) {
      tag.back() = c;
      int index = loop->find_tag(tag);
      if (index == -1)
        fail("Miller index tag not found: " + tag);
      indices.push_back(index);
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'H';
      col->label = alpha_up(c);
    }

    // M/ISYM and BATCH
    if (unmerged) {
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'Y';
      col->label = "M/ISYM";

      col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'B';
      col->label = "BATCH";
    }

    // other columns according to the spec
    bool column_added = false;
    for (const Entry& entry : spec_entries) {
      tag.replace(len, std::string::npos, entry.refln_tag);
      int index = loop->find_tag(tag);
      if (index == -1)
        continue;
      if (mtz.column_with_label(entry.col_label))
        continue;
      column_added = true;
      indices.push_back(index);
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = entry.dataset_id;
      col->type = entry.col_type;
      if (col->type == 's') {
        col->type = 'I';
        uses_status = true;
      }
      col->label = entry.col_label;
      if (verbose)
        out << "  " << tag << " -> " << col->label << '\n';
    }
    if (!column_added)
      fail(force_unmerged ? "Unmerged d" : "D", "ata not found in block ", rb.block.name);

    for (size_t i = 0; i != mtz.columns.size(); ++i) {
      mtz.columns[i].parent = &mtz;
      mtz.columns[i].idx = i;
    }
    mtz.nreflections = (int) loop->length();

    std::unique_ptr<UnmergedHklMover> hkl_mover;
    std::vector<std::pair<int,int>> batch_nums;
    if (unmerged) {
      hkl_mover.reset(new UnmergedHklMover(mtz.spacegroup));
      tag.replace(len, std::string::npos, "diffrn_id");
      int sweep_id_index = loop->find_tag(tag);
      tag.replace(len, std::string::npos, "pdbx_image_id");
      int image_id_index = loop->find_tag(tag);
      if (sweep_id_index == -1 || image_id_index == -1) {
        if (verbose)
          out << "No pdbx_image_id, setting BATCH to a dummy value.\n";
        auto batch = mtz.batches.emplace(mtz.batches.end());
        batch->number = 1;
        batch->set_dataset_id(1);
        batch->set_cell(mtz.cell);
        // FIXME should we set more properties in BATCH header?
      } else {
        if (verbose)
          out << "  " << tag << " & diffrn_id -> BATCH\n";
        // store sweep and frame numbers corresponding to reflections
        batch_nums.reserve(loop->length());
        for (size_t i = 0; i < loop->values.size(); i += loop->tags.size()) {
          int sweep_id = cif::as_int(loop->values[i + sweep_id_index], 0);
          const std::string& frame_str = loop->values[i + image_id_index];
          int frame = 1;
          if (!cif::is_null(frame_str))
            frame = (int) std::ceil(cif::as_number(frame_str));
          batch_nums.emplace_back(sweep_id, frame);
        }
        // store unique frame numbers
        std::map<int, std::set<int>> sets;
        for (const std::pair<int,int>& p : batch_nums)
          if (p.second >= 0)
            sets[p.first].insert(p.second);
        // add offset to frame numbers to make them unique
        std::map<int, int> offsets;
        int cap = 0;
        for (const auto& it : sets) {
          offsets.emplace(it.first, cap);
          cap += *--it.second.end() + 1100;
          cap -= cap % 1000;
        }
        for (std::pair<int,int>& p : batch_nums)
          if (p.first >= 0 && p.second >= 0)
            p.second += offsets.at(p.first);
        for (const auto& sweep_frames : sets) {
          Mtz::Batch batch;
          batch.set_dataset_id(sweep_frames.first);
          batch.set_cell(mtz.cell);
          int min_frame = *sweep_frames.second.begin();
          int max_frame = *--sweep_frames.second.end();
          int offset = offsets.at(sweep_frames.first);
          if (2 * sweep_frames.second.size() > size_t(max_frame - min_frame)) {
            // probably consecutive range, even if some frames are missing
            for (int n = min_frame; n <= max_frame; ++n) {
              batch.number = n + offset;
              mtz.batches.push_back(batch);
            }
          } else {
            for (int n : sweep_frames.second) {
              batch.number = n + offset;
              mtz.batches.push_back(batch);
            }
          }
        }
      }
    }

    // fill in the data
    mtz.data.resize(mtz.columns.size() * mtz.nreflections);
    size_t k = 0, row = 0;
    for (size_t i = 0; i < loop->values.size(); i += loop->tags.size()) {
      if (unmerged) {
        std::array<int, 3> hkl;
        for (int ii = 0; ii != 3; ++ii)
          hkl[ii] = cif::as_int(loop->values[i + indices[ii]]);
        int isym = hkl_mover->move_to_asu(hkl);
        for (int j = 0; j != 3; ++j)
          mtz.data[k++] = (float) hkl[j];
        mtz.data[k++] = (float) isym;
        mtz.data[k++] = batch_nums.empty() ? 1.f : (float) batch_nums[row++].second;
      } else {
        for (int j = 0; j != 3; ++j)
          mtz.data[k++] = (float) cif::as_int(loop->values[i + indices[j]]);
      }
      size_t j = 3;
      if (uses_status)
        mtz.data[k++] = status_to_freeflag(loop->values[i + indices[j++]]);
      for (; j != indices.size(); ++j) {
        const std::string& v = loop->values[i + indices[j]];
        if (cif::is_null(v)) {
          mtz.data[k] = (float) NAN;
        } else {
          mtz.data[k] = (float) cif::as_number(v);
          if (std::isnan(mtz.data[k]))
            out << "Value #" << i + indices[j] << " in the loop is not a number: "
                << v << '\n';
        }
        ++k;
      }
    }
    return mtz;
  }

private:
  static float status_to_freeflag(const std::string& str) {
    char c = str[0];
    if (c == '\'' || c == '"')
      c = str[1];
    if (c == 'o')
      return 1.f;
    if (c == 'f')
      return 0.f;
    return NAN;
  }
};

} // namespace gemmi
#endif
