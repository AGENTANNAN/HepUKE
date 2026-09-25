# ============================================================================
# 1605.03256v1  Observation of e+e- -> eta' J/psi
#   eta' -> eta pi+ pi- (eta -> gamma gamma)      and   eta' -> gamma pi+ pi-
#   J/psi -> l+ l-   (l = e, mu)
# Data: 14 c.m. energy points from sqrt(s) = 4.189 to 4.600 GeV (about 4.5 fb^-1
# of total integrated luminosity).
#
# BOSS-side spec: dataset preparation + event selection up to the final
# kinematic fit. The Born cross-section extraction (ISR correction, vacuum
# polarisation, upper limits) and the simultaneous mass fit are ROOT-level.
# ============================================================================

### ---------------------------------------------------------------------------
### Datasets — 14 energy points (sample keys chosen so that the tabulated
### luminosities match Table I of the paper; the narrow 4.189/4.208/4.217/
### 4.242/4.308/4.387/4.575 points belong to the 4260/4360 scans)
### ---------------------------------------------------------------------------
data_4189 = DatasetManager.real_data.find("703_4190scan")   # 4.189 GeV,   43.1 pb^-1
data_4208 = DatasetManager.real_data.find("703_4210scan")   # 4.208 GeV,   54.6 pb^-1
data_4217 = DatasetManager.real_data.find("703_4220scan")   # 4.217 GeV,   54.1 pb^-1
data_4226 = DatasetManager.real_data.find("703_4230")       # 4.226 GeV, 1047.3 pb^-1
data_4242 = DatasetManager.real_data.find("703_4245")       # 4.242 GeV,   55.6 pb^-1
data_4258 = DatasetManager.real_data.find("703_4260")       # 4.258 GeV,  825.7 pb^-1
data_4308 = DatasetManager.real_data.find("703_4310")       # 4.308 GeV,   44.9 pb^-1
data_4358 = DatasetManager.real_data.find("703_4360")       # 4.358 GeV,  539.8 pb^-1
data_4387 = DatasetManager.real_data.find("703_4390")       # 4.387 GeV,   55.2 pb^-1
data_4416 = DatasetManager.real_data.find("703_4420")       # 4.416 GeV, 1028.9 pb^-1
data_4467 = DatasetManager.real_data.find("703_4470")       # 4.467 GeV,  109.9 pb^-1
data_4527 = DatasetManager.real_data.find("703_4530")       # 4.527 GeV,  110.0 pb^-1
data_4575 = DatasetManager.real_data.find("703_4575")       # 4.575 GeV,   47.7 pb^-1
data_4600 = DatasetManager.real_data.find("703_4600")       # 4.600 GeV,  566.9 pb^-1

data_points = [data_4189, data_4208, data_4217, data_4226, data_4242, data_4258,
               data_4308, data_4358, data_4387, data_4416, data_4467, data_4527,
               data_4575, data_4600]

# Inclusive MC used for the background study; the paper generates it at
# sqrt(s) = 4.258, 4.416 and 4.600 GeV (Y(4260) decays, ISR production of the
# vector charmonium states, continuum hadrons and QED processes).
incMC_points = [DatasetManager.inclusive_mc.find("703_4260"),   # 4.258 GeV
                DatasetManager.inclusive_mc.find("703_4420"),   # 4.416 GeV
                DatasetManager.inclusive_mc.find("703_4600")]   # 4.600 GeV

### ---------------------------------------------------------------------------
### Decay cards — four independent reconstruction chains.
### The process is a continuum production of a final state with no intermediate
### charmonium, so the psi(4260) top mother (BESIII KKMC convention) is used.
### The Born cross section is fitted later as an incoherent sum of a
### Breit-Wigner for the psi(4160) and a polynomial continuum term.
### ---------------------------------------------------------------------------
# Chain I : eta' -> gamma pi+ pi-,  J/psi -> e+ e-   (1 photon, 4 tracks)
decay_card_chainI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Chain II : eta' -> gamma pi+ pi-,  J/psi -> mu+ mu-  (1 photon, 4 tracks)
decay_card_chainII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Chain III : eta' -> eta pi+ pi- (eta -> gamma gamma),  J/psi -> e+ e-
#             (2 photons, 4 tracks)
decay_card_chainIII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Chain IV : eta' -> eta pi+ pi- (eta -> gamma gamma),  J/psi -> mu+ mu-
#            (2 photons, 4 tracks)
decay_card_chainIV = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive signal MC — one sample per energy point for each of the four
### channels (the same card is used at every scan point, so 'Particle vpho' /
### a hard-coded energy is never written).
### ---------------------------------------------------------------------------
exMCs_chainI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etapJpsi_gamma2pi_ee"          # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_chainI
  config.cross_section = :default
end

exMCs_chainII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etapJpsi_gamma2pi_mumu"
  config.events        = 200_000
  config.decay_card    = decay_card_chainII
  config.cross_section = :default
end

exMCs_chainIII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etapJpsi_eta2pi_ee"
  config.events        = 200_000
  config.decay_card    = decay_card_chainIII
  config.cross_section = :default
end

exMCs_chainIV = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etapJpsi_eta2pi_mumu"
  config.events        = 200_000
  config.decay_card    = decay_card_chainIV
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — one Algorithm per chain (Rule T1: the four chains
### have different photon multiplicities and different kinematic-fit
### hypotheses). ECMS is NOT set via set_constant: this is a multi-energy
### cross-section scan, so the per-job c.m. energy is injected at run time.
### ---------------------------------------------------------------------------
# chain key -> (algorithm name, decay card, min photons, fit list).
# Chains III/IV additionally constrain M(gamma gamma) to the eta mass (5C fit),
# so need_eta is true for those two.
chains = [
  { name: "EtaPrimeJpsiGamma2PiEE",   card: decay_card_chainI,
    n_gam: ">=1", need_eta: false, fit: [:gamma, :pip, :pim, :lp, :lm] },
  { name: "EtaPrimeJpsiGamma2PiMuMu", card: decay_card_chainII,
    n_gam: ">=1", need_eta: false, fit: [:gamma, :pip, :pim, :lp, :lm] },
  { name: "EtaPrimeJpsiEta2PiEE",     card: decay_card_chainIII,
    n_gam: ">=2", need_eta: true,  fit: [:gamma, :gamma, :pip, :pim, :lp, :lm] },
  { name: "EtaPrimeJpsiEta2PiMuMu",   card: decay_card_chainIV,
    n_gam: ">=2", need_eta: true,  fit: [:gamma, :gamma, :pip, :pim, :lp, :lm] }
]

algorithms = {}
chains.each do |ch|
  alg = Algorithm.new(ch[:name])
  alg.set_header(["#{ch[:name]}Alg/#{ch[:name]}.h"])
     .set_alias({"std::vector<double>" => "Vdouble"})

  sel = Selection.new
  # Four good charged tracks with zero net charge
  sel.select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vz        10.0      # |R_z| < 10.0 cm along the beam direction
        Vr        1.0       # R_xy < 1.0 cm perpendicular to the beam
        nChrp     "==2"     # pi+/l+ pair
        nChrn     "==2"     # pi-/l- pair
        nNet      "==0"
      }
     # Good photons: at least 1 for eta' -> gamma pi+ pi-, at least 2 for
     # eta' -> eta pi+ pi- (eta -> gamma gamma)
     .select_photon {
        tdc_emc_start     0        # EMC cluster timing in [0, 700] ns
        tdc_emc_end       14
        angle_to_track    20.0     # > 20 deg from the nearest charged track
        energyThreshold_b 0.025    # > 25 MeV in the barrel (|cos(theta)| < 0.80)
        energyThreshold_e 0.050    # > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
        nGam              ch[:n_gam]
      }
     # Two pions (p < 0.8 GeV) and one lepton pair (p > 1.0 GeV) among the four
     # tracks; the e/mu separation is made by the E/p ratio (see notes).
     .pid(method: :probability) {
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :pion, against: [:kaon, :proton]
        npip "==1"
        npim "==1"
        nlp  "==1"                 # exactly one lepton in each charge
        nlm  "==1"
      }
     # Nominal kinematic fit to the hypothesis e+e- -> (gamma) pi+ pi- l+ l-.
     # Chain I/II : 4C on gamma pi+pi- l+l-
     # Chain III/IV: 5C, the additional constraint being the gamma gamma
     #               invariant mass fixed to the eta nominal mass
     .kinematic_fit(ch[:fit]) {
        nominal
        constrain_four_momentum
        if ch[:need_eta]
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        end
        chi2_cut 200               # loose BOSS cut; the paper's chi2_4C/chi2_5C < 40 is applied in ROOT
      }

  alg
    .note(:lepton_pid_epcut,
          "Charged-track classification: a track with momentum < 0.8 GeV is taken as a " \
          "pion candidate, a track with momentum > 1.0 GeV as a lepton candidate. " \
          "Electron/muon separation uses the E/p ratio of the EMC energy deposit and the " \
          "MDC momentum: E/p > 0.8 for electrons, E/p < 0.4 for muons. Both tracks of the " \
          "lepton pair must be of the same species.")
    .note(:best_photon_combination,
          "For eta' -> gamma pi+ pi- the photon candidate giving the minimum chi2_4C is " \
          "retained when the event contains more than one good photon; for " \
          "eta' -> eta pi+ pi- the photon pair giving the minimum chi2_5C is chosen.")
    .note(:background_veto,
          "ISR psi(3686) background and FSR from the leptons are removed by requiring " \
          "M(pi+ pi- J/psi) and the recoil mass M_recoil(pi+ pi-) to lie outside the " \
          "regions 3.65-3.71 GeV (3.67-3.71 for the eta pi+ pi- channel) and " \
          "3.05-3.15 GeV (3.65-3.69) respectively. e+e- -> eta psi(3686) -> eta pi+ pi- J/psi " \
          "and e+e- -> pi+ pi- psi(3686) -> pi+ pi- eta J/psi are the reactions vetoed.")
    .note(:signal_region,
          "The eta' is selected in the window 0.94 < M(gamma(eta) pi+ pi-) < 0.98 GeV/c^2, " \
          "with sidebands 0.90-0.94 and 0.98-1.02 GeV/c^2; the J/psi is selected in " \
          "3.07 < M(l+ l-) < 3.13 GeV/c^2. The invariant-mass distributions are fitted " \
          "simultaneously over the four channels in ROOT (MC signal shape convolved with " \
          "a Gaussian plus a linear background), and the statistical significance and the " \
          "90% C.L. upper limits are extracted there.")
    .note(:radiative_correction,
          "The Born cross section sigma_B = N_obs / (L_int (1+delta) |1+Pi|^2 sum_i eps_i B_i) " \
          "uses the ISR radiative correction factor (1+delta) from the radiator function and " \
          "the vacuum-polarisation factor |1+Pi|^2; both are computed outside the BOSS selection.")
    .with_decay_card(ch[:card])
    .apply(sel)

  algorithms[ch[:name]] = alg
end

### ---------------------------------------------------------------------------
### Execution — real data + inclusive MC + the four exclusive signal samples
### ---------------------------------------------------------------------------
root_files = {}
root_files["chainI"]   = algorithms["EtaPrimeJpsiGamma2PiEE"].execute_on(
                             data_points + incMC_points + exMCs_chainI)
root_files["chainII"]  = algorithms["EtaPrimeJpsiGamma2PiMuMu"].execute_on(
                             data_points + incMC_points + exMCs_chainII)
root_files["chainIII"] = algorithms["EtaPrimeJpsiEta2PiEE"].execute_on(
                             data_points + incMC_points + exMCs_chainIII)
root_files["chainIV"]  = algorithms["EtaPrimeJpsiEta2PiMuMu"].execute_on(
                             data_points + incMC_points + exMCs_chainIV)
