### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")        # psi(3686) data, 2.7B events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # psi(3686) inclusive MC

# Decay card for psi(3686) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(2S)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  e+  e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card for psi(3686) -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(2S)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  mu+  mu-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples
exMC_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_3686_pipijpsi_ee"
    config.related_dataset = psip_data
    config.events = 100000
    config.decay_card = decay_card_ee
    config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_3686_pipijpsi_mumu"
    config.related_dataset = psip_data
    config.events = 100000
    config.decay_card = decay_card_mumu
    config.cross_section = :default
end

### Event selection (BOSS) — Channel I: J/psi -> e+ e- ###
alg_ee = Algorithm.new("PsippToPiPiJpsiEE")
alg_ee.set_header(["PsippToPiPiJpsiEEAlg/PsippToPiPiJpsiEE.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:helix_correction,
          "helix-parameter track correction applied to all four final-state particles before the 4C kinematic fit; \
           systematic uncertainty estimated from efficiency difference between with/without correction")
        .note(:tracking_efficiency_correction,
          "2D tracking efficiency correction in (p_T, cos_theta) for pi+ pi- derived from control sample \
           psi(3686) -> pi+ pi- J/psi with J/psi -> e+ e-")
        .note(:cp_odd_observable,
          "CP-odd observable q3 = p_e+ · (p_pi+ - p_pi-) × p_e+ · (p_pi+ × p_pi-) / |p_pi+ × p_pi-| \
           computed in e+e- CM frame; events with |q3| < 0.24 (GeV/c)^3 kept; \
           efficiency reweighting N_obs^i / epsilon_i applied per q3 bin; \
           A_CP = (N(q3>0)-N(q3<0)) / (N(q3>0)+N(q3<0))")
        .note(:continuum_suppression,
          "continuum background (e+e- -> pi+pi- J/psi at sqrt(s)=3.65 GeV) suppressed with \
           |cos_theta_{pi,l}| > 0.95 and m(pi+pi-) < 0.32 GeV/c^2; \
           continuum events normalized by luminosity ratio × (s_3650/s_3686) for background subtraction")
        .note(:muon_energy_threshold,
          "paper identifies muon when EMC deposited energy < 0.45 GeV; \
           DSL identifies muon when EMC <= treat_as_electron_if_energy_above (1.0 GeV here); \
           gap region [0.45, 1.0] GeV may cause minor discrepancy")

sel_ee = Selection.new
sel_ee.select_track{
            cos_theta   0.93       # |cos(theta)| < 0.93
            Vz          10.0       # |Vz| < 10 cm in beam direction
            Vr          1.0        # Vr < 1 cm in transverse plane
            nChrp       "==2"      # Exactly 2 positive tracks
            nChrn       "==2"      # Exactly 2 negative tracks
            nNet        "==0"      # Net charge zero
        }
        .pid(method: :probability) do
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 1.0
            identify :pion, against: [:kaon, :proton]
            npip   "==1"     # One pi+
            npim   "==1"     # One pi-
            nep    "==1"     # One e+
            nem    "==1"     # One e-
        end
        .kinematic_fit([:pip, :pim, :ep, :em]) do
            nominal
            vertex_fit([0, 1, 2, 3])   # primary vertex fit on all four tracks
            constrain_four_momentum    # 4C energy-momentum constraint
            invariant_mass_of(:ep, :em).within(3.087, 3.107)  # |m(e+e-) - M_J/psi| < 0.01 GeV/c^2
            chi2_cut 200               # loose chi2 cut; tight cut in ROOT
        end

alg_ee.with_decay_card(decay_card_ee).apply(sel_ee)

### Event selection (BOSS) — Channel II: J/psi -> mu+ mu- ###
alg_mumu = Algorithm.new("PsippToPiPiJpsiMuMu")
alg_mumu.set_header(["PsippToPiPiJpsiMuMuAlg/PsippToPiPiJpsiMuMu.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:helix_correction,
            "helix-parameter track correction applied to all four final-state particles before the 4C kinematic fit; \
             systematic uncertainty estimated from efficiency difference between with/without correction")
          .note(:tracking_efficiency_correction,
            "2D tracking efficiency correction in (p_T, cos_theta) for pi+ pi- derived from control sample \
             psi(3686) -> pi+ pi- J/psi with J/psi -> e+ e-; muon tracking efficiency referenced from pions")
          .note(:cp_odd_observable,
            "CP-odd observable q3 = p_e+ · (p_pi+ - p_pi-) × p_e+ · (p_pi+ × p_pi-) / |p_pi+ × p_pi-| \
             computed in e+e- CM frame; events with |q3| < 0.24 (GeV/c)^3 kept; \
             efficiency reweighting N_obs^i / epsilon_i applied per q3 bin; \
             A_CP = (N(q3>0)-N(q3<0)) / (N(q3>0)+N(q3<0))")
          .note(:continuum_suppression,
            "continuum background (e+e- -> pi+pi- J/psi at sqrt(s)=3.65 GeV) suppressed with \
             |cos_theta_{pi,l}| > 0.95 and m(pi+pi-) < 0.32 GeV/c^2; \
             continuum events normalized by luminosity ratio × (s_3650/s_3686) for background subtraction")
          .note(:muon_energy_threshold,
            "paper identifies muon when EMC deposited energy < 0.45 GeV; \
             DSL identifies muon when EMC <= treat_as_electron_if_energy_above (1.0 GeV here); \
             gap region [0.45, 1.0] GeV may cause minor discrepancy")

sel_mumu = Selection.new
sel_mumu.select_track{
            cos_theta   0.93       # |cos(theta)| < 0.93
            Vz          10.0       # |Vz| < 10 cm in beam direction
            Vr          1.0        # Vr < 1 cm in transverse plane
            nChrp       "==2"      # Exactly 2 positive tracks
            nChrn       "==2"      # Exactly 2 negative tracks
            nNet        "==0"      # Net charge zero
        }
        .pid(method: :probability) do
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 1.0
            identify :pion, against: [:kaon, :proton]
            npip   "==1"     # One pi+
            npim   "==1"     # One pi-
            nmup   "==1"     # One mu+
            nmum   "==1"     # One mu-
        end
        .kinematic_fit([:pip, :pim, :mup, :mum]) do
            nominal
            vertex_fit([0, 1, 2, 3])   # primary vertex fit on all four tracks
            constrain_four_momentum    # 4C energy-momentum constraint
            invariant_mass_of(:mup, :mum).within(3.087, 3.107)  # |m(mu+mu-) - M_J/psi| < 0.01 GeV/c^2
            chi2_cut 200               # loose chi2 cut; tight cut in ROOT
        end

alg_mumu.with_decay_card(decay_card_mumu).apply(sel_mumu)

# Execute algorithms on datasets
root_files_ee = alg_ee.execute_on([psip_data, psip_incMC, exMC_ee])
root_files_mumu = alg_mumu.execute_on([psip_data, psip_incMC, exMC_mumu])