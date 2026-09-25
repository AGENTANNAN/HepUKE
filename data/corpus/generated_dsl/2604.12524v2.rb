# ============================================================================
# ψ(2S) → γ χ_c1, χ_c1 → π+π− η'   @ 3.686 GeV   (BOSS part: dataset prep +
# event selection up to and including the kinematic fit)
#   Channel 1 : η' → γ π+π−            (Dalitz)
#   Channel 2 : η' → π+π− η,  η → γγ
# Two independent final states / selection chains → two Algorithm objects.
# ============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(2S) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # matching inclusive MC

# --- Decay card, channel 1: η' → γ π+π− (Dalitz) ---
decay_card_ch1 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.000 pi+ pi- eta' PHSP;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Decay card, channel 2: η' → π+π− η with η → γγ ---
decay_card_ch2 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.000 pi+ pi- eta' PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive MC, 500k events for each η' decay mode ---
exMC_ch1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_etap_gammapipi"   # η' → γ π+π−
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_ch1
  config.cross_section   = :default
end

exMC_ch2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_etap_pipieta"     # η' → π+π− η, η → γγ
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_ch2
  config.cross_section   = :default
end

# ============================================================================
### Event selection — Channel 1 : η' → γ π+π− ###
# ============================================================================
alg_name_1 = "ChiC1EtaPGammaPiPi"
alg_ch1 = Algorithm.new(alg_name_1)
alg_ch1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_ch1 = Selection.new
sel_ch1.select_track {                       # charged-track quality cuts
          cos_theta 0.93                     # |cosθ| < 0.93
          Vz        10.0                     # |Vz| < 10 cm
          Vr        10.0                     # Vr < 10 cm
          nChrp     "==2"                    # exactly two positive tracks
          nChrn     "==2"                    # exactly two negative tracks
          nNet      "==0"                    # net charge zero
        }
        .select_photon {                     # photon selection
          tdc_emc_start     0                # EMC timing window [0, 14] (×700 ns)
          tdc_emc_end       14
          angle_to_track    10.0             # opening angle to nearest charged track > 10°
          energyThreshold_b 0.100            # E > 100 MeV (barrel)
          energyThreshold_e 0.100            # E > 100 MeV (endcap)
          nGam              ">=2"            # channel 1: at least two photons
        }
        # no PID is applied — all charged tracks are treated as pions
        .assign({:chrgp => :pip, :chrgn => :pim})
        # nominal 4C kinematic fit to γ γ π+ π− π+ π−
        # (loose χ² cut here; tight χ²_4C < 40 applied later in ROOT)
        .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) {
          nominal
          constrain_four_momentum
          invariant_mass_of(:gamma, :pip, :pim).within(0.9128, 1.0028)          # |M(γπ+π−) − m_η'| < 45 MeV/c²
          invariant_mass_of(:gamma, :pip, :pim, :pip, :pim).larger_than(3.4)    # M(γπ+π−π+π−) > 3.4 GeV/c²
          chi2_cut 200
        }

alg_ch1.note(:combination_selection,
             "require exactly one γ π+ π− combination with |M(γπ+π−) − m_η'| < 45 MeV/c² " \
             "and M(γπ+π−π+π−) > 3.4 GeV/c²; combination counting is not expressible in DSL")

alg_ch1.with_decay_card(decay_card_ch1).apply(sel_ch1)
alg_ch1.execute_on([psip_data, psip_incMC, exMC_ch1])

# ============================================================================
### Event selection — Channel 2 : η' → π+π− η, η → γγ ###
# ============================================================================
alg_name_2 = "ChiC1EtaPPiPiEta"
alg_ch2 = Algorithm.new(alg_name_2)
alg_ch2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_ch2 = Selection.new
sel_ch2.select_track {                       # charged-track quality cuts
          cos_theta 0.93
          Vz        10.0
          Vr        10.0
          nChrp     "==2"
          nChrn     "==2"
          nNet      "==0"
        }
        .select_photon {                     # photon selection
          tdc_emc_start     0                # EMC timing window [0, 14] (×700 ns)
          tdc_emc_end       14
          angle_to_track    10.0             # opening angle to nearest charged track > 10°
          energyThreshold_b 0.100            # E > 100 MeV (barrel)
          energyThreshold_e 0.100            # E > 100 MeV (endcap)
          nGam              ">=3"            # channel 2: at least three photons
        }
        # no PID is applied — all charged tracks are treated as pions
        .assign({:chrgp => :pip, :chrgn => :pim})
        # 4C kinematic fit to γ γ γ π+ π− π+ π−
        .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
          constrain_four_momentum
          chi2_cut 200
        }
        # nominal 5C fit: 4C + η mass constraint on a γγ pair
        # (loose χ² cut; tight χ²_5C < 40 applied later in ROOT)
        .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
          nominal
          constrain_four_momentum
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # 5th constraint: M(γγ) = m_η
          invariant_mass_of(:pip, :pim, :gamma, :gamma).within(0.9368, 0.9788)   # |M(π+π−γγ) − m_η'| < 21 MeV/c²
          invariant_mass_of(:gamma, :gamma, :pip, :pim, :pip, :pim).larger_than(3.4) # M(γγπ+π−π+π−) > 3.4 GeV/c²
          chi2_cut 200
        }

alg_ch2.note(:combination_selection,
             "require exactly one π+ π− γ γ combination with |M(π+π−γγ) − m_η'| < 21 MeV/c² " \
             "and M(γγπ+π−π+π−) > 3.4 GeV/c²; combination counting is not expressible in DSL")

alg_ch2.with_decay_card(decay_card_ch2).apply(sel_ch2)
alg_ch2.execute_on([psip_data, psip_incMC, exMC_ch2])