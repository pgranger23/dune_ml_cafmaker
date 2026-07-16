# Design Proposal: SPINE ML Reconstruction Integration in DUNE CAFs

Following the feedback on the initial schema draft, this proposal defines a clean, modular, and algorithm-agnostic standard record structure for the Common Analysis Files (CAFs). 

---

## 1. Summary of Feedback & Design Solutions

### 1. Rename `dlp` to `spine`
* **Action:** Rename all code occurrences of `dlp`/`DLP` to `spine`/`SPINE` (specifically variables, branches, and class names).
* **Rationale:** SPINE is the name of the Deep Learning/ML reconstruction stack used in the analysis chain.

### 2. Move Detector-Specific Metrics to `SRLArRecoParticle` (Renamed from `SRPFP`)
* **Action:** Remove TPC drift-specific and cathode-specific metrics (`is_cathode_crosser`, `is_time_contained`, `cathode_offset`) from the generic `SRRecoParticle` class in `common`.
* **Action:** Merge `SRSpineInfo` into `SRLArRecoParticle` and rename the container vector `pfps` inside `SRFDInt` to `lar_particles`.
* **Rationale:** Keeping detector-specific or stack-specific reconstruction metrics inside `common` makes the generic particle representation cluttered. Placing them under `SRLArRecoParticle` puts them cleanly under the detector-specific `fd.hd`/`fd.vd` branches of the CAF (e.g. `rec.fd.hd.pandora.lar_particles`).

### 3. Simplify `SRInteraction::part` (Direct Vector of Particles)
* **Action:** Completely eliminate the `SRRecoParticlesBranch` wrapper class.
* **Action:** Change `SRInteraction::part` to be a direct `std::vector<SRRecoParticle>`.
* **Rationale:** Because interactions are already split by reconstruction stack at the top-level (e.g. `rec.common.ixn.spine` vs `rec.common.ixn.pandora`), all particles inside an interaction belong to the same reconstruction stack. Therefore, having a wrapper class containing separate vectors is redundant. A direct vector of `SRRecoParticle` is the cleanest design.

### 4. Group PID and Energy Estimators (`SREnergy` under `SRPIDProbability`)
* **Action:** Move the hypothesis-specific energy estimators (`calo`, `csda`, `mcs` kinetic energy) into the PID class.
* **Action:** Define an `SRPIDHypothesis` class that groups the PDG ID, the classification score, the $\chi^2$ value, and the corresponding `SREnergy` estimates.
* **Action:** Avoid `std::map<int, SREnergy>` due to a known bug in `gen_srproxy` where the parser crashes when traversing `typedef` internals of STL map templates. Instead, store them as a `std::vector<SRPIDHypothesis>` and provide map-like helper methods (`PIDScore(int pdg)`, `Chi2(int pdg)`, `Energy(int pdg)`) for a clean analyzer interface.
* **Action:** Keep measured, hypothesis-independent values (`depositions`, `start_dedx`) directly in `SRRecoParticle`.

### 5. Remove `SRClusteringInfo` & Merge with `SRLArRecoParticle`
* **Action:** Remove the `SRClusteringInfo` class entirely from standard record headers.
* **Action:** Add `size` and `num_fragments` directly to `SRLArRecoParticle`.
* **Rationale:** Clustering and fragment metrics are highly dependent on the hit/voxelization layout and reconstruction stack, which are detector/reconstruction stack specific features.

---

## 2. Detailed Class Layouts

### A. PID & Energy Layouts (`duneanaobj/StandardRecord/SRPIDProbability.h`)

```cpp
namespace caf
{
  class SREnergy
  {
    public:
      static constexpr float NaN = std::numeric_limits<float>::signaling_NaN();

      float calo = NaN;                  ///< Calorimetric kinetic energy [GeV]
      float csda = NaN;                  ///< CSDA range-based kinetic energy [GeV]
      float mcs = NaN;                   ///< MCS deflection-based kinetic energy [GeV]
  };

  class SRPIDHypothesis
  {
    public:
      static constexpr float NaN = std::numeric_limits<float>::signaling_NaN();

      int      pdg = 0;                  ///< PDG code for the hypothesis (e.g. 13 for muon)
      float    pid_score = NaN;          ///< Probability score for this species [0 - 1]
      float    chi2 = NaN;               ///< Fit chi2 value under this hypothesis
      SREnergy energy;                   ///< Estimated energy under this hypothesis
  };

  class SRPIDProbability
  {
    public:
      static constexpr float NaN = std::numeric_limits<float>::signaling_NaN();

      std::vector<SRPIDHypothesis>  hyp;              ///< Hypothesis list (PDG, score, chi2, and energy)
      std::vector<float>            primary_scores;   ///< Probability scores for primary vs. secondary (or track-like classifications)

      // Map-like helper methods
      float PIDScore(int pdg) const;
      float Chi2(int pdg) const;
      const SREnergy& Energy(int pdg) const;
  };
}
```

### B. Updated Interaction Member (`duneanaobj/StandardRecord/SRInteraction.h`)

```cpp
namespace caf
{
  class SRInteraction
  {
    public:
      // ... (existing members like id, vtx, dir, nuhyp, Enu, truth, contained)

      /// Collections of reconstructed particles associated with this interaction
      std::vector<SRRecoParticle> part;
  };
}
```

### C. Generic Particle Extensions (`duneanaobj/StandardRecord/SRRecoParticle.h`)

```cpp
namespace caf
{
  class SRRecoParticle
  {
    public:
      // ... (existing members like p, start, end, contained, score, truth)

      // Generic kinematic extensions
      SRVector3D  start_dir;              ///< Reconstructed direction vector at start
      SRVector3D  end_dir;                ///< Reconstructed direction vector at end
      float       length = NaN;           ///< Reconstructed path length [cm]

      // Measured energy quantities (independent of particle species hypothesis)
      float       depositions = NaN;      ///< Total raw energy deposition [GeV]
      float       start_dedx = NaN;       ///< dE/dx near the start of the path [MeV/cm]

      // Energy estimates under the assigned (best) hypothesis
      float       calo_ke = NaN;          ///< Assigned calorimetric kinetic energy [GeV]
      float       csda_ke = NaN;          ///< Assigned CSDA range-based kinetic energy [GeV]
      float       mcs_ke = NaN;           ///< Assigned MCS deflection-based kinetic energy [GeV]

      // Modular structures
      SRPIDProbability pid_prob;          ///< Detailed PID probabilities and hypothesis energy estimates
  };
}
```

### D. Detector-Specific Particle Flow Object (`duneanaobj/StandardRecord/SRLArRecoParticle.h`)

```cpp
namespace caf
{
  class SRLArRecoParticle : public SRRecoObjBase
  {
    public:
      static constexpr float NaN = std::numeric_limits<float>::signaling_NaN();

      // Hit counts per view
      int nhits_U = 0;                    ///< Number of hits in the U wire plane
      int nhits_V = 0;                    ///< Number of hits in the V wire plane
      int nhits_W = 0;                    ///< Number of hits in the W wire plane
      int nhits_3D = 0;                   ///< Number of 3D hits/space points

      // LArPfoHierarchyFeatureTool
      float daughter_parent_hit_ratio = NaN; ///< Ratio of hits in daughters vs. parent
      float ndaughters_hit_3d = NaN;      ///< Number of 3D hits in all daughters

      // LArThreeDChargeFeatureTool
      float charge_end_fraction = NaN;    ///< Fraction of charge deposited near the end of the track
      float charge_fractional_spread = NaN; ///< Fractional spread of the charge

      // LArThreeDLinearFitFeatureTool
      float diff_straight_line_mean = NaN; ///< Mean difference from a straight line fit
      float line_length = NaN;            ///< Length of the reconstructed line fit
      float max_fit_gap_length = NaN;     ///< Maximum gap between hits in the fit
      float sliding_linear_fit_rms = NaN; ///< RMS of the sliding linear fit

      // LArThreeDOpeningAngleFeatureTool
      float angle_diff_3d = NaN;          ///< Opening angle/direction difference in 3D

      // LArThreeDPCAFeatureTool
      float secondary_pca_ratio = NaN;    ///< Ratio of secondary to primary principal component axis
      float tertiary_pca_ratio = NaN;     ///< Ratio of tertiary to primary principal component axis

      // LArThreeDVertexDistanceFeatureTool
      float vertex_distance = NaN;        ///< Distance to the reconstructed vertex

      // Track/Shower BDT score
      float track_score = NaN;            ///< BDT classifier score separating tracks from showers [0 - 1]

      // TPC and drift-specific reco metrics (Moved from generic particle level)
      bool  is_cathode_crosser = false;   ///< Does the particle cross the cathode plane?
      bool  is_time_contained = false;    ///< Is the particle time-contained within the drift window?
      float cathode_offset = NaN;         ///< Distance/offset to the cathode [cm]

      // Clustering and fragment details (Merged from SRClusteringInfo)
      int64_t size = -1;                  ///< Number of hits/voxels/points in the cluster
      int64_t num_fragments = -1;         ///< Number of sub-clusters or fragments

      // SPINE-specific reconstruction details (Merged from SRSpineInfo)
      bool        is_valid = false;             ///< Is the SPINE reconstruction valid for this particle
      LArTPCPid   pid = kLArTPCUnknown;         ///< SPINE particle species classification
      LArTPCShape shape = kLArTPCUnknownShape;  ///< SPINE particle shape classification
      float       start_straightness = NaN;     ///< Straightness of the track near start
      float       axial_spread = NaN;           ///< Axial voxel spread
      float       directional_spread = NaN;     ///< Directional voxel spread
  };
}
```

---

## 3. Map-Like Helper Methods Implementation

To make it easy for analyzers to access hypothesis-specific metrics by PDG without exposing ROOT to the template dictionary compilation crashes, we implement helper methods in `SRPIDProbability.cxx`:

```cpp
#include "duneanaobj/StandardRecord/SRPIDProbability.h"

namespace caf
{
  float SRPIDProbability::PIDScore(int pdg) const
  {
      for (const auto& h : hyp) {
          if (h.pdg == pdg) return h.pid_score;
      }
      return NaN;
  }

  float SRPIDProbability::Chi2(int pdg) const
  {
      for (const auto& h : hyp) {
          if (h.pdg == pdg) return h.chi2;
      }
      return NaN;
  }

  const SREnergy& SRPIDProbability::Energy(int pdg) const
  {
      for (const auto& h : hyp) {
          if (h.pdg == pdg) return h.energy;
      }
      // Return a static default SREnergy containing NaN values if not found
      static const SREnergy default_energy;
      return default_energy;
  }
}
```
