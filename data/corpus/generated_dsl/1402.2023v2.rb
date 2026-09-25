# ============================================================================
#  ψ(3686) → γ χ_cJ (J = 1, 2),  χ_cJ → η' K+ K-
#    Mode I  : η' → γ ρ0,  ρ0 → π+ π-          → final state γγ K+ K- π+ π-
#    Mode II : η' → η π+ π-,  η → γ γ           → final state γγγ K+ K- π+ π-
#  (BOSS part: dataset preparation + event selection up to the nominal 4C fit)
# ============================================================================

### ----------------------------- Datasets --------------------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Decay card — Mode I : ψ(2S) → γ χ_c1, χ_c1 → η' K+ K-, η' → γ ρ0, ρ0 → π+ π-
decay_card_modeI = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.000 eta' K+ K- VSS;
  Enddecay

  Decay eta'
  1.000 gamma rho0 PHSP;
  Enddecay

  Decay rho0
  1.000 pi+ pi- VSS;
  Enddecay

  End
DECAYCARD

# Decay card — Mode II : ψ(2S) → γ χ_c1, χ_c1 → η' K+ K-, η' → η π+ π-, η → γ γ
decay_card_modeII = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.000 eta' K+ K- VSS;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive MC for each mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_chic1_etapkk_modeI"
  config.related_dataset = psip_data          # associated real dataset
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_chic1_etapkk_modeII"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### -------------------- Mode I : η' → γ ρ0, ρ0 → π+ π- -------------------- ###
alg_name_I = "EtapGamRhoKK"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV

sel_modeI = Selection.new
  .select_track {
    cos_theta  0.93     # |cosθ| < 0.93
    Vz         10.0     # |Vz| < 10 cm
    Vr         1.0      # Vr < 1 cm
    nChrp      "==2"    # exactly two positive tracks
    nChrn      "==2"    # exactly two negative tracks
    nNet       "==0"    # net charge zero
  }
  .select_photon {
    angle_to_track     5.0     # photon more than 5° from any charged track
    energyThreshold_b  0.025   # barrel energy threshold 25 MeV
    energyThreshold_e  0.050   # endcap energy threshold 50 MeV
    nGam               ">=2"   # at least two photons (radiative γ + η' daughter γ)
  }
  .pid(method: :probability) {
    prob_cut  0.001
    identify  :kaon, against: [:pion, :proton]   # identify K+ and K-
    nkp       ">=1"                              # at least one K+
    nkm       ">=1"                              # at least one K-
  }
  .assign({:chrgp => :pip, :chrgn => :pim})      # remaining tracks assigned as π+ / π-
  .remove([:kp <= :pip, :km <= :pim])            # drop overlapping K/π assignments
  # Nominal 4C kinematic fit on γγ K+ K- π+ π-
  .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut  40
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

### -------------- Mode II : η' → η π+ π-, η → γ γ -------------------------- ###
alg_name_II = "EtapEtaPiPiKK"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV

sel_modeII = Selection.new
  .select_track {
    cos_theta  0.93     # |cosθ| < 0.93
    Vz         10.0     # |Vz| < 10 cm
    Vr         1.0      # Vr < 1 cm
    nChrp      "==2"    # exactly two positive tracks
    nChrn      "==2"    # exactly two negative tracks
    nNet       "==0"    # net charge zero
  }
  .select_photon {
    angle_to_track     5.0     # photon more than 5° from any charged track
    energyThreshold_b  0.025   # barrel energy threshold 25 MeV
    energyThreshold_e  0.050   # endcap energy threshold 50 MeV
    nGam               ">=3"   # at least three photons (radiative γ + η → γγ)
  }
  .pid(method: :probability) {
    prob_cut  0.001
    identify  :kaon, against: [:pion, :proton]   # identify K+ and K-
    nkp       ">=1"                              # at least one K+
    nkm       ">=1"                              # at least one K-
  }
  .assign({:chrgp => :pip, :chrgn => :pim})      # remaining tracks assigned as π+ / π-
  .remove([:kp <= :pip, :km <= :pim])            # drop overlapping K/π assignments
  # Nominal 4C kinematic fit on γγγ K+ K- π+ π-
  .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut  50
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])