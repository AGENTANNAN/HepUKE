# =====================================================================
# Inclusive J/psi yield measurement at BESIII
#   N_J/psi = (N_sel - N_bg) / ( eps_trig * eps_data^{psi(3686)} * f_cor )
# Inclusive counting analysis: no exclusive final state, no kinematic fit.
# The same inclusive selection chain is applied to the J/psi data, psi(3686)
# data, the inclusive MC and both exclusive MC samples.
# =====================================================================

### ----------------------------- Datasets ----------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi data (2009 + 2012) @ 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC @ 3.097 GeV
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC @ 3.686 GeV

### ---------------------------- Decay cards --------------------------- ###
# (1) psi(3686) -> pi+ pi- J/psi , with J/psi -> anything
decay_card_psip_to_pipi_jpsi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    End
DECAYCARD

# (2) J/psi -> anything , produced directly at rest at the J/psi energy
decay_card_jpsi_inclusive = <<~DECAYCARD
    Decay J/psi
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    End
DECAYCARD

### -------------------------- Exclusive MC ----------------------------- ###
# 1M-event exclusive MC samples
exMC_psip_to_pipi_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_to_pipi_jpsi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip_to_pipi_jpsi
  config.cross_section   = :default
end

exMC_jpsi_inclusive = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_inclusive"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi_inclusive
  config.cross_section   = :default
end

### ---------------- Inclusive selection (ECMS = 3.097 GeV) ---------------- ###
alg_incl_name = "InclusiveJpsi"
alg_incl = Algorithm.new(alg_incl_name)
alg_incl.set_header(["#{alg_incl_name}Alg/#{alg_incl_name}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        # Event-level cuts without a dedicated selection keyword -> recorded as notes
        .note(:track_momentum_cut, "each good charged track required to have p < 2.0 GeV/c")
        .note(:visible_energy,     "event visible energy E_vis > 1.0 GeV")
        .note(:bhabha_dimuon_veto, "when exactly two charged tracks are present, require both track momenta < 1.5 GeV/c (Bhabha / dimuon veto)")
        .note(:emc_deposit_cut,    "per-track EMC energy deposit < 1 GeV")

sel_incl = Selection.new
sel_incl
  .select_track {                        # inclusive charged-track selection
      cos_theta 0.93                     # |cos(theta)| < 0.93
      Vr        1.0                      # Vr < 1 cm (transverse plane)
      Vz        15.0                     # |Vz| < 15 cm
      nTot      ">=2"                    # at least two charged tracks
  }
  .select_photon {                       # inclusive photon selection
      tdc_emc_start     0                # EMC cluster timing 0 < T <= 700 ns
      tdc_emc_end       14
      energyThreshold_b 0.025            # E > 25 MeV in the barrel (|cos(theta)| < 0.83)
      energyThreshold_e 0.050            # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.93)
  }

# Same inclusive chain applied to J/psi and psi(3686) data, inclusive MC and both exclusive MC
alg_incl.with_decay_card(decay_card_jpsi_inclusive).apply(sel_incl)
alg_incl.execute_on([jpsi_data, jpsi_incMC, psip_data, psip_incMC,
                     exMC_psip_to_pipi_jpsi, exMC_jpsi_inclusive])

### ------------ psi(3686) tag side : pi+ pi- recoil at 3.686 GeV ------------ ###
alg_tag_name = "PsipTag"
alg_tag = Algorithm.new(alg_tag_name)
alg_tag.set_header(["#{alg_tag_name}Alg/#{alg_tag_name}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_tag = Selection.new
sel_tag
  .select_track {                        # tag: >= 2 oppositely charged pions
      cos_theta 0.93
      Vr        1.0
      Vz        15.0
      nChrp     ">=1"
      nChrn     ">=1"
  }
  .assign({:chrgp => :pip, :chrgn => :pim})
  .remove(:pip) { condition "three_momentum_of(:pip) > 0.4" }  # keep p < 0.4 GeV/c
  .remove(:pim) { condition "three_momentum_of(:pim) > 0.4" }  # no further cut on remaining tracks / showers

alg_tag.with_decay_card(decay_card_psip_to_pipi_jpsi).apply(sel_tag)
alg_tag.execute_on([psip_data, exMC_psip_to_pipi_jpsi])