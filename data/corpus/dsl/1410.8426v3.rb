# =============================================================================
# BESIII arXiv:1410.8426v3
# Search for the weak decays J/psi -> D_s^(*)- e+ nu_e + c.c.
#
# Data: 2.25 x 10^8 J/psi events at sqrt(s) = 3.097 GeV (BOSS sample 708_3097).
#       A 478 pb^-1 sample at sqrt(s) = 4.009 GeV (703_4009) is used only as a
#       control sample, psi(4040) -> D_s^+ D_s^-, to determine the systematic
#       uncertainty of the D_s reconstruction efficiency; it is reconstructed
#       with a separate hadronic single-tag and is NOT processed by the signal
#       selection below.
#
# The D_s^- is reconstructed in four decay modes:
#   (1) K+ K- pi-        (2) K+ K- pi- pi0
#   (3) K_S0 K-          (4) K_S0 K+ pi- pi-
# and the D_s*- candidate from D_s*- -> D_s- gamma.  Together with the positron
# and the undetected neutrino this completes J/psi -> D_s^(*)- e+ nu_e.
# Charge conjugation is implied throughout.
#
# Two signal channels (D_s and D_s*) x four D_s decay modes -> eight Algorithm
# objects (Rule T1: independent decay modes get separate Algorithm objects).
# Scope: dataset preparation + event selection (BOSS side) only.
# =============================================================================

### ------------------------------- Datasets -------------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 2.25e8 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC

# Control sample for the D_s reconstruction efficiency systematic uncertainty
# (psi(4040) -> D_s^+ D_s^- at 4.009 GeV); not processed by the signal selection.
psi4040_data  = DatasetManager.real_data.find("703_4009")
psi4040_incMC = DatasetManager.inclusive_mc.find("703_4009")

### ---------------------------- Decay cards -------------------------------- ###
# Per-mode description of the D_s^- decay and of the event topology.
#   :body   -- D_s^- daughters written into the decay card
#   :sub    -- extra sub-decay block (K_S0 -> pi+ pi- or pi0 -> gamma gamma)
#   :ngam   -- minimum photon multiplicity in the D_s and D_s* channels
#   :eneut  -- maximum total energy (GeV) of the extra neutral particles
mode_table = [
  { tag: "KKPi", body: "K+ K- pi-", sub: nil,
    ngam: { ds: nil, dsstar: 1 }, eneut: { ds: 0.20, dsstar: 0.30 } },

  { tag: "KKPiPi0", body: "K+ K- pi- pi0",
    sub: "Decay pi0\n1.0000 gamma gamma PHSP;\nEnddecay",
    ngam: { ds: 2, dsstar: 3 }, eneut: { ds: 0.15, dsstar: 0.20 } },

  { tag: "KsK", body: "K_S0 K-",
    sub: "Decay K_S0\n1.0000 pi+ pi- PHSP;\nEnddecay",
    ngam: { ds: nil, dsstar: 1 }, eneut: { ds: 0.20, dsstar: 0.30 } },

  { tag: "KsKPiPi", body: "K_S0 K+ pi- pi-",
    sub: "Decay K_S0\n1.0000 pi+ pi- PHSP;\nEnddecay",
    ngam: { ds: nil, dsstar: 1 }, eneut: { ds: 0.20, dsstar: 0.30 } }
]

channels = { ds: "Ds", dsstar: "DsStar" }

# Build the EvtGen decay card for one (channel, D_s decay mode) combination.
# The J/psi -> D_s^(*) e nu decay is generated with the dedicated weak-interaction
# generator of the paper, i.e. assuming the c -> s charged-current process.
build_card = lambda do |channel, body, sub|
  lines = ["Decay J/psi"]
  lines << (channel == :dsstar ? "1.0000 D_s*- e+ nu_e PHSP;" : "1.0000 D_s- e+ nu_e PHSP;")
  lines << "Enddecay"
  lines << ""
  if channel == :dsstar
    lines << "Decay D_s*-"
    lines << "1.0000 gamma D_s- PHSP;"
    lines << "Enddecay"
    lines << ""
  end
  lines << "Decay D_s-"
  lines << "1.0000 #{body} PHSP;"
  lines << "Enddecay"
  lines << ""
  if sub
    lines << sub
    lines << ""
  end
  lines << "End"
  lines.join("\n") + "\n"
end

### ------------------- Exclusive signal MC (100k per mode) ------------------ ###
# The paper generates 100,000 events for each D_s^- decay mode and each channel.
exMC = {}
mode_table.each do |mode|
  channels.each_key do |channel|
    exMC[[channel, mode[:tag]]] = DatasetManager.create_exclusive_mc do |config|
      config.sample_name     = "JpsiTo#{channels[channel]}ENu_#{mode[:tag]}_exclusive_mc"
      config.related_dataset = jpsi_data
      config.events          = 100_000
      config.decay_card      = build_card.call(channel, mode[:body], mode[:sub])
      config.cross_section   = :default
    end
  end
end

### --------------------- Event selection (BOSS) ---------------------------- ###
# One Algorithm per (channel, D_s decay mode) pair.  The neutrino is undetected,
# so the decay chain is completed by partial_miss instead of a 4C kinematic fit.
# RecID layout (level-order traversal of the decay card), independent of the mode:
#   0 J/psi | 1 D_s^(*)- | 2 e+ | 3 nu_e | 4 gamma (only in the D_s* channel) | 5 D_s- | ...
# so the single missed particle is recID 3 in both channels.
jobs = []

mode_table.each do |mode|
  channels.each_key do |channel|
    alg_name = "JpsiTo#{channels[channel]}ENu#{mode[:tag]}"
    alg = Algorithm.new(alg_name)
    alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)

    sel = Selection.new

    # Charged tracks: |cos(theta)| < 0.93, within +-10 cm of the IP along the beam
    # and within +-1 cm transverse to the beam.  The exact track multiplicities
    # implement the requirement that no charged particle other than those from the
    # D_s^- and the positron candidate is present in the event.  For the
    # K_S0 K+ pi- pi- mode the final state holds six charged tracks (3 positive,
    # 3 negative), for the other three modes four (2 positive, 2 negative).
    if mode[:tag] == "KsKPiPi"
      sel.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==3"
        nChrn     "==3"
        nTot      "==6"
      end
    else
      sel.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nTot      "==4"
      end
    end

    # Photons: E > 25 MeV in the barrel (|cos(theta)| < 0.80), E > 50 MeV in the
    # endcap (0.86 < |cos(theta)| < 0.92), more than 20 degrees from any charged
    # track, and EMC timing 0 < T < 700 ns.  The required multiplicity counts the
    # pi0 daughters and, in the D_s* channel, the transition photon.
    sel.select_photon do
      tdc_emc_start     0
      tdc_emc_end       14        # 14 x 50 ns = 700 ns
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    20.0
      ngam = mode[:ngam][channel]
      nGam ">=#{ngam}" if ngam
    end

    # PID from the combined dE/dx and TOF probabilities:
    #   e+-  P(e) > 0.001, P(e) > P(K), P(e) > P(pi)
    # The positron candidate is identified first; it is then removed from the
    # charged-track lists so that it cannot be reused as a pion or kaon.
    sel.pid(method: :probability) do
      prob_cut 0.001
      identify :electron, against: [:pion, :kaon]   # the single positron candidate
      nep "==1"
    end

    # 0.80 < E/p < 1.05 for the positron candidate, where E is the energy deposited
    # in the EMC and p the momentum measured in the MDC.
    sel.for_each(:ep) do
      define(:e_over_p) { eraw / p }
      where { (e_over_p < 0.80) | (e_over_p > 1.05) }
      remove
    end
    sel.remove([:ep <= :chrgp])

    # The remaining charged tracks are identified as kaons and pions:
    #   K+-  P(K) > 0.001, P(K) > P(pi)
    #   pi+- P(pi) > 0.001, P(pi) > P(K)
    # The K_S0 daughters are treated by the secondary vertex fit below.
    sel.pid(method: :probability) do
      prob_cut 0.001
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon]
      case mode[:tag]
      when "KKPi", "KKPiPi0"
        nkp  "==1"
        nkm  "==1"
        npim "==1"
      when "KsK"
        nkm  "==1"
        npip "==1"
        npim "==1"
      when "KsKPiPi"
        nkp  "==1"
        npip "==1"
        npim "==3"
      end
    end

    # K_S0 -> pi+ pi-: a primary and a secondary vertex fit are performed, the
    # decay length must exceed twice its fit error, and the invariant mass of the
    # pair must satisfy 0.487 < M(pi+ pi-) < 0.511 GeV/c^2.  Multiple K_S0
    # candidates are allowed in one event.
    if mode[:tag].start_with?("Ks")
      sel.secondary_vertex_fit([:pip, :pim]) do
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      end
    end

    # pi0 -> gamma gamma: kinematic fit constraining the photon-pair invariant mass
    # to the nominal pi0 mass, retaining the minimum-chi2 combination with
    # chi2 < 100 and 0.115 < M(gamma gamma) < 0.150 GeV/c^2.  Candidates with both
    # photons in the endcap are rejected because of the poor resolution there.
    if mode[:tag] == "KKPiPi0"
      sel.kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 100
        npi0 ">=1"
      end
    end

    # The neutrino is undetected; every other particle of the chain is
    # reconstructed.  The missing four-momentum and U_miss = E_miss - |p_miss|
    # are formed from the recoil four-momentum and are the signal variables.
    sel.partial_miss([3]) do
      # recID 3 = nu_e (massless, undetected); recID 0 is the J/psi top mother.
    end

    alg.note(:generator_model,
             "The signal J/psi -> D_s^(*)- e+ nu_e is generated with a dedicated " \
             "generator assuming the process is dominated by the weak interaction " \
             "(the c -> s charged-current process); hadronization and quark-spin-flip " \
             "effects are ignored.")
       .note(:electron_pid,
             "Positron identification requires P(e) > 0.001, P(e) > P(K) and " \
             "P(e) > P(pi) from the combined dE/dx and TOF probabilities, plus the " \
             "energy-momentum ratio 0.80 < E/p < 1.05, where E is the energy deposited " \
             "in the EMC and p the momentum measured by the MDC. The E/p condition is " \
             "applied as an explicit per-track filter on the positron candidate list.")
       .note(:pion_kaon_pid,
             "Kaon candidates require P(K) > 0.001 and P(K) > P(pi); pion candidates " \
             "require P(pi) > 0.001 and P(pi) > P(K).")
       .note(:ks0_selection,
             "K_S0 candidates are built from all pairs of oppositely charged tracks " \
             "assumed to be pions without any PID requirement, with the " \
             "interaction-point requirements relaxed to 20 cm along the beam " \
             "direction. A primary and a secondary vertex fit are performed and the " \
             "K_S0 decay length must be more than two times larger than its fit error. " \
             "The invariant mass, computed from the track parameters after the " \
             "secondary vertex fit, must satisfy 0.487 < M(pi+ pi-) < " \
             "0.511 GeV/c^2. Multiple K_S0 candidates are allowed in one event.")
       .note(:pi0_selection,
             "pi0 candidates are reconstructed from photon pairs; the combination " \
             "with the minimum chi^2 of the mass-constrained fit that satisfies " \
             "chi^2 < 100 and 0.115 < M(gamma gamma) < 0.150 GeV/c^2 is retained. " \
             "Candidates with both photons in the endcap regions are excluded.")
       .note(:extra_neutral_energy,
             "To reduce background from misidentified events with extra photons, the " \
             "total energy of the extra neutral particles is required to be less than " \
             "#{mode[:eneut][channel]} GeV for the #{mode[:tag]} mode of the " \
             "#{channel == :ds ? 'D_s-' : 'D_s*-'} channel. The criterion was chosen by " \
             "optimising S/sqrt(B), with S from the signal MC and B from the inclusive " \
             "MC, and is applied as an event-level sum over the unused photon " \
             "candidates.")
       .note(:charged_veto,
             "The signal candidate must contain a positron track, and events that " \
             "include charged particles other than those from the D_s^- and the " \
             "positron candidate are vetoed. This is implemented by the exact " \
             "nChrp/nChrn/nTot conditions combined with the PID multiplicity " \
             "conditions.")
       .note(:ds_mass_window,
             "The invariant mass of the D_s^- candidate is required to lie within " \
             "three times the mode-specific mass resolution (+-3 sigma around the " \
             "central value); the window differs between the four decay modes because " \
             "the resolutions differ, and the numeric windows are applied in the " \
             "offline analysis.")
       .note(:missing_momentum,
             "|p_miss| > 50 MeV is required, with " \
             "p_miss = p_J/psi - p_D_s^(*)- - p_e+; this suppresses the background " \
             "from J/psi hadronic decays in which a pion is misidentified as a " \
             "positron.")
       .note(:signal_variable,
             "The events are extracted with U_miss = E_miss - |p_miss|, where " \
             "E_miss = E_J/psi - E_D_s^(*)- - E_e+. For correctly identified " \
             "semileptonic decays U_miss peaks at zero. A simultaneous unbinned " \
             "maximum-likelihood fit over the four D_s decay modes is performed in " \
             "-0.2 < U_miss < 0.2 GeV/c^2, with the signal described by a Gaussian plus " \
             "a Crystal Ball function and the background by the inclusive MC shape; it " \
             "is carried out offline on the NTuple.")
       .note(:control_sample,
             "A 478 pb^-1 data sample collected at sqrt(s) = 4.009 GeV " \
             "(psi(4040) -> D_s^+ D_s^-) is used to study the systematic uncertainty of " \
             "the D_s reconstruction efficiency: one D_s is tagged in eight hadronic " \
             "decay modes and the other D_s is reconstructed in the same way as in this " \
             "analysis. The differences between the D_s reconstruction efficiencies in " \
             "MC and data are quoted as the systematic uncertainty. This study uses its " \
             "own tag-side reconstruction and is not part of the signal selection.")
       .note(:upper_limit,
             "No significant excess of signal is observed in either channel. Bayesian " \
             "upper limits at the 90% C.L. are set on the number of signal events " \
             "(244 for J/psi -> D_s^- e+ nu_e and 335 for J/psi -> D_s*^- e+ nu_e), " \
             "corresponding to B(J/psi -> D_s^- e+ nu_e + c.c.) < 1.3 x 10^-6 and " \
             "B(J/psi -> D_s*^- e+ nu_e + c.c.) < 1.8 x 10^-6. This is a statistical " \
             "analysis performed offline.")

    if channel == :dsstar
      alg.note(:dsstar_reconstruction,
               "D_s*- candidates are reconstructed by combining the D_s^- candidate " \
               "with an additional photon candidate and requiring the invariant-mass " \
               "difference m(D_s^- gamma) - m(D_s^-) to lie in " \
               "0.125 < Delta M < 0.150 GeV/c^2. No requirement is made to select the " \
               "best D_s^- or D_s*- candidate and multiple candidates are allowed; " \
               "multiple-candidate events occur in about 0.2% of the events per mode " \
               "in MC.")
    else
      alg.note(:multi_candidate,
               "No requirement is made to select the best D_s^- candidate and multiple " \
               "D_s^- candidates are allowed in one event; multiple-candidate events " \
               "occur in about 0.1% of the events per mode in MC.")
    end

    alg.with_decay_card(build_card.call(channel, mode[:body], mode[:sub])).apply(sel)

    jobs << [alg, exMC[[channel, mode[:tag]]]]
  end
end

### ------------------------------ Execution --------------------------------- ###
# Each algorithm is run on the J/psi data, the matching inclusive MC and its own
# exclusive signal MC.  The 4.009 GeV control sample is not processed here.
jobs.each { |alg, mc| alg.execute_on([jpsi_data, jpsi_incMC, mc]) }
