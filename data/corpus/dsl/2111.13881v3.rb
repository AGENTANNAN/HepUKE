DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) -> pi+ pi- J/psi, J/psi -> e+ e- e+ e-
decay_card_4e = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi JPIPI;
    Enddecay
    Decay J/psi
    1.000 e+ e- e+ e- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card: psi(3686) -> pi+ pi- J/psi, J/psi -> e+ e- mu+ mu-
decay_card_2e2mu = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi JPIPI;
    Enddecay
    Decay J/psi
    1.000 e+ e- mu+ mu- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card: psi(3686) -> pi+ pi- J/psi, J/psi -> mu+ mu- mu+ mu-
decay_card_4mu = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi JPIPI;
    Enddecay
    Decay J/psi
    1.000 mu+ mu- mu+ mu- PHSP;
    Enddecay
    End
DECAYCARD

exMC_4e    = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_4e";     c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_4e;    c.cross_section = :default }
exMC_2e2mu = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_2e2mu";  c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_2e2mu; c.cross_section = :default }
exMC_4mu   = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_4mu";    c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_4mu;   c.cross_section = :default }

# Mode I: J/psi -> e+ e- e+ e-
alg_4e = Algorithm.new("JpsiTo4e")
alg_4e.set_header(["JpsiTo4eAlg/JpsiTo4e.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .note(:soft_pion_cut,
         "Two oppositely charged soft pions with p < 0.45 GeV/c selected for pi+ pi- J/psi tag. " \
         "Combination with M_rec(pi+ pi-) closest to J/psi mass within 5-sigma window retained. " \
         "sigma = 4.0 +/- 0.2 MeV/c^2 for this channel.")
       .note(:photon_conversion_veto,
         "Photon conversion finder used: R_xy is the distance from IP to the conversion point " \
         "of the lower-momentum e+ e- pair. Delta_p = (p_e_h - p_e_l) > 1.0 GeV/c required " \
         "to separate signal from gamma-conversion background.")
       .note(:helix_correction,
         "Helix parameter correction applied to simulated tracks before kinematic fit, " \
         "derived from psi(3686) -> pi+ pi- J/psi, J/psi -> l+ l- control sample.")

sel_4e = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        ">=6"
    nNet        "==0"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    prob_cut 0.001
  end
  .kinematic_fit([:pip, :pim, :ep, :em, :ep, :em]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_4e.with_decay_card(decay_card_4e).apply(sel_4e)

# Mode II: J/psi -> e+ e- mu+ mu-
alg_2e2mu = Algorithm.new("JpsiTo2e2mu")
alg_2e2mu.set_header(["JpsiTo2e2muAlg/JpsiTo2e2mu.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .note(:soft_pion_cut,
            "Two oppositely charged soft pions with p < 0.45 GeV/c selected for pi+ pi- J/psi tag. " \
            "Combination with M_rec(pi+ pi-) closest to J/psi mass within 5-sigma window retained. " \
            "sigma = 2.8 +/- 0.5 MeV/c^2 for this channel.")
          .note(:muon_pid_cuts,
            "Muon candidates: prob(mu) > prob(e) AND prob(mu) > prob(K). " \
            "EMC deposited energy in [0.1, 0.3] GeV. " \
            "Electron candidates: prob(e) > prob(mu) AND prob(e) > prob(K).")
          .note(:photon_conversion_veto,
            "Photon conversion finder with R_xy and Delta_p > 1.0 GeV/c cut, same as 4e channel.")
          .note(:competing_bkg_veto,
            "Events with chi2_4C(J/psi -> e+ e- mu+ mu-) < chi2_4C(J/psi -> e+ e- pi+ pi-) retained " \
            "to suppress J/psi -> e+ e- pi+ pi- background. ROOT-level cut on competing chi2.")
          .note(:helix_correction,
            "Helix parameter correction applied to simulated tracks before kinematic fit.")

sel_2e2mu = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        ">=6"
    nNet        "==0"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    prob_cut 0.001
  end
  # Nominal signal fit: e+ e- mu+ mu-
  .kinematic_fit([:pip, :pim, :ep, :em, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end
  # Competing hypothesis: e+ e- pi+ pi- (Rule T2)
  .kinematic_fit([:pip, :pim, :ep, :em, :pip, :pim]) do
    constrain_four_momentum
  end

alg_2e2mu.with_decay_card(decay_card_2e2mu).apply(sel_2e2mu)

# Mode III: J/psi -> mu+ mu- mu+ mu-
alg_4mu = Algorithm.new("JpsiTo4mu")
alg_4mu.set_header(["JpsiTo4muAlg/JpsiTo4mu.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .note(:soft_pion_cut,
         "Two oppositely charged soft pions with p < 0.45 GeV/c selected for pi+ pi- J/psi tag. " \
         "Combination with M_rec(pi+ pi-) closest to J/psi mass retained.")
       .note(:muon_pid_cuts,
         "Muon candidates: prob(mu) > prob(e) AND prob(mu) > prob(K). " \
         "EMC deposited energy in [0.1, 0.3] GeV. Additional muon counter depth requirements applied.")
       .note(:competing_bkg_veto,
         "Events with chi2_4C(J/psi -> mu+ mu- mu+ mu-) < chi2_4C(J/psi -> pi+ pi- pi+ pi-) retained " \
         "to suppress J/psi -> pi+ pi- pi+ pi- background. ROOT-level cut on competing chi2.")
       .note(:helix_correction,
         "Helix parameter correction applied to simulated tracks before kinematic fit.")
       .note(:peaking_background,
         "Two peaking background events from J/psi -> pi+ pi- pi+ pi- observed in inclusive MC; " \
         "accounted for in the fit at ROOT level.")

sel_4mu = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        ">=6"
    nNet        "==0"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    prob_cut 0.001
  end
  # Nominal signal fit: mu+ mu- mu+ mu-
  .kinematic_fit([:pip, :pim, :lp, :lm, :lp, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end
  # Competing hypothesis: pi+ pi- pi+ pi- (Rule T2)
  .kinematic_fit([:pip, :pim, :pip, :pim, :pip, :pim]) do
    constrain_four_momentum
  end

alg_4mu.with_decay_card(decay_card_4mu).apply(sel_4mu)

datasets = [psip_data, psip_incMC, exMC_4e, exMC_2e2mu, exMC_4mu]
alg_4e.execute_on(datasets)
alg_2e2mu.execute_on(datasets)
alg_4mu.execute_on(datasets)