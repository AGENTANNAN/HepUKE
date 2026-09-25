### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Signal decay card: ψ(2S) → π+π− J/ψ, J/ψ → K+K−
decay_card_KK = <<~DECAYCARD
  Decay psi(2S)
  1.000  pi+  pi-  J/psi    PHSP;
  Enddecay

  Decay J/psi
  1.000  K+  K-    VSS;
  Enddecay

  End
DECAYCARD

# Reference decay card: ψ(2S) → π+π− J/ψ, J/ψ → μ+μ−
decay_card_mumu = <<~DECAYCARD
  Decay psi(2S)
  1.000  pi+  pi-  J/psi    PHSP;
  Enddecay

  Decay J/psi
  1.000  mu+  mu-    PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Exclusive MC: 750k events for each J/ψ decay mode
exMC_KK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pipijpsi_KK"
  config.related_dataset = psip_data
  config.events          = 750_000
  config.decay_card      = decay_card_KK
  config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pipijpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 750_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common selection for both channels: four charged tracks, no photon selection
# (all-charged final state), no PID at selection level.
common_selection = Selection.new
  .select_track {
    cos_theta 0.80     # |cosθ| < 0.80
    Vz        10.0     # |Vz| < 10 cm
    Vr        1.0      # Vr < 1 cm
    nChrp     "==2"    # exactly two positively charged tracks
    nChrn     "==2"    # exactly two negatively charged tracks
    nNet      "==0"    # net charge 0
  }
  # No photon selection applied (all-charged final state).
  .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: all positives → π+, all negatives → π−
  .kinematic_fit([:pip, :pip, :pim, :pim]) {  # common vertex fit on all four charged tracks
    nominal                                   # flagged nominal (endpoint of the BOSS selection)
    vertex_fit([0, 1, 2, 3])                  # no explicit nC constraint; vertex only
    chi2_cut 200                              # χ² < 200
  }

# --- Signal channel: J/ψ → K+K− ---
alg_KK = Algorithm.new("PsipPipiJpsiKK")
alg_KK.set_header(["PsipPipiJpsiKKAlg/PsipPipiJpsiKK.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .note(:momentum_based_kpi_separation,
            "no PID at selection level; K/π separation performed by momentum: charged tracks with p < 1.0 GeV/c treated as pions, p > 1.2 GeV/c treated as kaons")
      .note(:emc_total_energy_window,
            "total EMC energy window 0.3-2.5 GeV imposed as a later BOSS-level selection (no dedicated DSL method)")
      .note(:electron_veto,
            "electron veto E_dep/p < 0.8 imposed on charged tracks as a later BOSS-level selection (no dedicated DSL method)")

sel_KK = common_selection.dup
alg_KK.with_decay_card(decay_card_KK).apply(sel_KK)

# --- Reference channel: J/ψ → μ+μ− ---
alg_mumu = Algorithm.new("PsipPipiJpsiMumu")
alg_mumu.set_header(["PsipPipiJpsiMumuAlg/PsipPipiJpsiMumu.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .note(:momentum_based_kpi_separation,
              "no PID at selection level; K/π separation performed by momentum: charged tracks with p < 1.0 GeV/c treated as pions, p > 1.2 GeV/c treated as kaons")
        .note(:emc_total_energy_window,
              "total EMC energy window 0.3-2.5 GeV imposed as a later BOSS-level selection (no dedicated DSL method)")
        .note(:electron_veto,
              "electron veto E_dep/p < 0.8 imposed on charged tracks as a later BOSS-level selection (no dedicated DSL method)")

sel_mumu = common_selection.dup
alg_mumu.with_decay_card(decay_card_mumu).apply(sel_mumu)

# --- Execute on real data, inclusive MC and the channel-specific exclusive MC ---
root_files_KK   = alg_KK.execute_on([psip_data, psip_incMC, exMC_KK])
root_files_mumu = alg_mumu.execute_on([psip_data, psip_incMC, exMC_mumu])