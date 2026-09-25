# =============================================================================
# BOSS DSL:  e+e- -> D*+D*-  and  e+e- -> D*+D-
# at 28 c.m. energy points spanning 4.085 - 4.600 GeV.
# Only the D*+ -> pi+ D0 -> pi+ K- pi+  side is reconstructed; the D*-/D- side
# is inferred from energy-momentum conservation (via the 4C kinematic fit).
# =============================================================================

### ---------------------------- Datasets ---------------------------------- ###
# The 28 c.m. energy points (4.085 - 4.600 GeV), named [BOSS]_[Ecms(MeV)]
energy_points = [
  "703_4090", "705_4130", "705_4160", "703_4190", "703_4200",
  "703_4210", "703_4220", "703_4230", "703_4237", "703_4245",
  "703_4246", "703_4260", "703_4270", "703_4280", "705_4290",
  "703_4310", "705_4315", "705_4340", "703_4360", "705_4380",
  "703_4390", "705_4400", "703_4420", "705_4440", "703_4470",
  "703_4530", "703_4575", "703_4600"
]

data_points  = energy_points.map { |name| DatasetManager.real_data.find(name) }     # matching real data
incMC_points = energy_points.map { |name| DatasetManager.inclusive_mc.find(name) }   # matching inclusive MC

### --------------------------- Decay cards -------------------------------- ###
# e+e- -> D*+ D*-  (both D* are open-charm; the D*+ side is reconstructed)
decay_card_dst_dstbar = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D*-  PHSP;
    Enddecay

    Decay D*+
    1.000 pi+ D0  PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+  PHSP;
    Enddecay

    Decay D*-
    1.000 pi- anti-D0  PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi-  PHSP;
    Enddecay

    End
DECAYCARD

# e+e- -> D*+ D-  (the D- side is not reconstructed)
decay_card_dst_d = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D-  PHSP;
    Enddecay

    Decay D*+
    1.000 pi+ D0  PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+  PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi-  PHSP;
    Enddecay

    End
DECAYCARD

### -------------------- Exclusive MC (100k events, per mode) -------------- ###
# One 100k-event signal MC per energy point for each of the two modes.
exMC_dst_dstbar = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dst_dstbar"
  config.events        = 100_000
  config.decay_card    = decay_card_dst_dstbar
  config.cross_section = :default
end

exMC_dst_d = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dst_d"
  config.events        = 100_000
  config.decay_card    = decay_card_dst_d
  config.cross_section = :default
end

### --------------- Shared event selection (both signal modes) ------------- ###
selection = Selection.new
    .select_track {
      cos_theta  0.93        # |cos(theta)| < 0.93
      Vr         1.0         # Vr < 1 cm (transverse plane)
      Vz         10.0        # |Vz| < 10 cm (beam axis)
      nChrp      ">=2"       # at least two positive charged tracks
      nChrn      ">=1"       # at least one negative charged track
    }
    # no photon requirement -> select_photon is intentionally not called
    .pid(method: :probability) {
      prob_cut 0.001                                             # probability cut = 0.001
      identify :kaon, against: [:pion, :proton]                  # K+ and K-
      identify :pion, against: [:kaon, :proton]                  # pi+ and pi-
      nkm  ">=1"                                                 # at least one K-
      npip ">=2"                                                 # at least two pi+
    }
    .remove(:km) { condition "three_momentum_of(:km) < 0.3" }    # K- candidates must have p > 0.3 GeV/c
    # split the pi+ candidates by momentum:
    #   pi+_L (soft, from D*+ -> pi+ D0)   : p < 0.3 GeV/c
    #   pi+_H (hard, from D0 -> K- pi+)    : p > 0.3 GeV/c
    .assign({:pip => :pip_L})
    .remove(:pip_L) { condition "three_momentum_of(:pip_L) > 0.3" }
    .assign({:pip => :pip_H})
    .remove(:pip_H) { condition "three_momentum_of(:pip_H) < 0.3" }
    # nominal kinematic fit on the K- pi+ pi+ system
    .kinematic_fit([:km, :pip_H, :pip_L]) {
      nominal                                                            # nominal fit (corrected 4-momenta saved)
      vertex_fit([0, 1])                                                 # constrain K- and pi+_H to a common (D0) vertex
      constrain_four_momentum                                            # 4C energy-momentum constraint
      invariant_mass_of(:km, :pip_H).within(1.8454, 1.8852)             # M(K- pi+_H) within +-3sigma of the D0 mass
      invariant_mass_of(:km, :pip_H).constrain_to_nominal_mass_of(:D0)  # constrain M(K- pi+_H) to the nominal D0 mass
      chi2_cut 200                                                       # chi2 < 200; smallest-chi2 combination chosen automatically
    }

### ---------------- Algorithms (one per signal mode) ---------------------- ###
# Mode I: e+e- -> D*+ D*-
alg_dst_dstbar = Algorithm.new("DstDstbar")
alg_dst_dstbar.set_header(["DstDstbarAlg/DstDstbar.h"])
              .set_constant({ "ECMS" => [:double, 4.26] })
              .note(:momentum_split, "The two pi+ are separated by momentum: pi+_L (from D*+ -> pi+ D0) with p < 0.3 GeV/c and pi+_H (from D0 -> K- pi+) with p > 0.3 GeV/c; at least one of each is required, and the D0 mass constraint is built from the K- pi+_H pair.")
              .note(:ecms_per_point, "Analysis runs over 28 c.m. energies (4.085-4.600 GeV); the 4C kinematic fit uses the per-run CMS energy, so the ECMS constant here is only a nominal fallback.")
alg_dst_dstbar.with_decay_card(decay_card_dst_dstbar).apply(selection)

# Mode II: e+e- -> D*+ D-
alg_dst_d = Algorithm.new("DstD")
alg_dst_d.set_header(["DstDAlg/DstD.h"])
          .set_constant({ "ECMS" => [:double, 4.26] })
          .note(:momentum_split, "The two pi+ are separated by momentum: pi+_L (from D*+ -> pi+ D0) with p < 0.3 GeV/c and pi+_H (from D0 -> K- pi+) with p > 0.3 GeV/c; at least one of each is required, and the D0 mass constraint is built from the K- pi+_H pair.")
          .note(:ecms_per_point, "Analysis runs over 28 c.m. energies (4.085-4.600 GeV); the 4C kinematic fit uses the per-run CMS energy, so the ECMS constant here is only a nominal fallback.")
alg_dst_d.with_decay_card(decay_card_dst_d).apply(selection.dup)

### --------------------------- Execute ------------------------------------ ###
root_files_dst_dstbar = alg_dst_dstbar.execute_on(data_points + incMC_points + exMC_dst_dstbar)
root_files_dst_d      = alg_dst_d.execute_on(data_points + incMC_points + exMC_dst_d)