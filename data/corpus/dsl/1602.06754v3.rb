# BESIII: psi -> Xi- Xi+ and Sigma(1385)-+ Sigma(1385)+- (psi = J/psi, psi(3686))
# Single baryon tag method: the Xi-(Sigma(1385)-) is fully reconstructed in
# pi-+ Lambda with Lambda -> p pi-, and the antibaryon partner is inferred from
# the recoil mass of the pi-+ Lambda system. Data: (223.7 +- 1.4)e6 J/psi and
# (106.4 +- 0.9)e6 psi(3686) events.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) data sample
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(3686) inclusive MC

# Signal decay cards. One million MC events are generated for each mode; the Xi
# and Sigma(1385) decays in the signal modes are simulated inclusively according
# to the PDG branching fractions, and the events are generated with the measured
# angular distribution parameter alpha.
decay_card_jpsi_xi = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi- anti-Xi+                        PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+                    PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                         HypWK;
  Enddecay

  End
DECAYCARD

decay_card_jpsi_sigma = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma- anti-Sigma+                  PHSP;
  Enddecay

  Decay Sigma-
  1.0000 Lambda0 pi-                         PHSP;
  Enddecay

  Decay anti-Sigma+
  1.0000 anti-Lambda0 pi+                    PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-                              HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                         HypWK;
  Enddecay

  End
DECAYCARD

decay_card_psip_xi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Xi- anti-Xi+                        PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+                    PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                         HypWK;
  Enddecay

  End
DECAYCARD

decay_card_psip_sigma = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma- anti-Sigma+                  PHSP;
  Enddecay

  Decay Sigma-
  1.0000 Lambda0 pi-                         PHSP;
  Enddecay

  Decay anti-Sigma+
  1.0000 anti-Lambda0 pi+                    PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi-                              HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+                         HypWK;
  Enddecay

  End
DECAYCARD

exMC_jpsi_xi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_to_Xim_Xip"
  config.related_dataset = jpsi_data
  config.events          = 1000000
  config.decay_card      = decay_card_jpsi_xi
  config.cross_section   = :default
end

exMC_jpsi_sigma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_to_Sigma1385m_Sigma1385p"
  config.related_dataset = jpsi_data
  config.events          = 1000000
  config.decay_card      = decay_card_jpsi_sigma
  config.cross_section   = :default
end

exMC_psip_xi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_Xim_Xip"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_psip_xi
  config.cross_section   = :default
end

exMC_psip_sigma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_Sigma1385m_Sigma1385p"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_psip_sigma
  config.cross_section   = :default
end

# Common selection template: single baryon tag via pi-+ Lambda with
# Lambda -> p pi-; the antibaryon is inferred from the recoil mass of the
# pi-+ Lambda system (partial reconstruction, no kinematic fit).
def build_tag_selection
  sel = Selection.new
  sel.select_track {
        cos_theta 0.93      # tracks within the MDC angular coverage
        Vz        10.0      # |Vz| < 10 cm (good helix fit at the IP)
        Vr        1.0        # |Vr| < 1 cm
        nChrp     ">=1"     # Xi- Xi+ : >=1 positive track (p from Lambda)
        nChrn     ">=2"     # Xi- Xi+ : >=2 negative tracks (p-bar and pi-)
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]     # pi-+ from Xi/Sigma(1385) decay
        identify :proton, against: [:kaon, :pion]     # p from Lambda decay
        npip ">=2"        # at least two charged pions
        nprp ">=1"        # at least one proton
      }
     # The proton and one pion of the tag are consumed by the Lambda
     .remove([:prp <= :chrgp])
     .remove([:pim <= :chrgn])
     # Lambda -> p pi- : vertex fit to all p pi- combinations (chi2 < 500),
     # |M(p pi-) - M_Lambda| < 6 MeV/c^2
     .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
  sel
end

########################################################################
# Algorithm I: J/psi -> Xi- Xi+ (single baryon tag: Xi- -> pi- Lambda)
########################################################################
alg_jpsi_xi = Algorithm.new("JpsiToXimXip")
alg_jpsi_xi.set_header(["JpsiToXimXipAlg/JpsiToXimXip.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })
           .with_decay_card(decay_card_jpsi_xi)
           .note(:baryon_tag,
                 "Single baryon tag: only the Xi- is reconstructed, as pi- Lambda " \
                 "with Lambda -> p pi-, and the anti-Xi+ partner is extracted from " \
                 "the recoil mass M_recoil(pi- Lambda) = sqrt((E_CM - E_piLambda)^2 - " \
                 "p_piLambda^2). The charge-conjugate tag mode Xi+ -> pi+ anti-Lambda " \
                 "is covered by the conjugate multiplicity requirement " \
                 "(nChrp >= 2, nChrn >= 1).")
           .note(:tag_mass_window,
                 "The pi-+ Lambda system is required to lie in the Xi- signal region " \
                 "[1.312, 1.332] GeV/c^2 (J/psi) / [1.308, 1.338] GeV/c^2 " \
                 "(psi(3686)); the candidate with the minimum " \
                 "|M(pi Lambda) - M_Xi(nominal)| is selected.")
           .note(:anti_baryon_signal_region,
                 "Signal events are those in the M_recoil(pi-+ Lambda) peak around the " \
                 "anti-Xi+ nominal mass; the recoil-mass spectrum is fitted with the " \
                 "simulated MC signal shape convolved with a Gaussian for the " \
                 "data/MC resolution difference and a second-order polynomial for the " \
                 "combinatorial background.")
           .note(:psip_to_pipim_jpsi_veto,
                 "For psi(3686) -> Xi- Xi+ an additional requirement " \
                 "|M_recoil(pi+ pi-) - M_J/psi| > 0.005 GeV/c^2 suppresses the " \
                 "background psi(3686) -> pi+ pi- J/psi (applied downstream in ROOT).")
           .note(:continuum_background,
                 "Continuum e+e- -> Xi- Xi+ background is estimated from data at " \
                 "3.08 GeV (300 nb^-1) and 3.65 GeV (44 pb^-1) with the same " \
                 "selection; no obvious peaking structure survives and the " \
                 "contribution is negligible.")
           .note(:angular_distribution,
                 "dN/d(cos theta) ~ 1 + alpha cos^2 theta, with theta the angle between " \
                 "the baryon and the e+ beam in the CM frame; fitted downstream in 16 " \
                 "cos theta bins over [-0.8, 0.8].")
           .note(:background_processes,
                 "Potential peaking backgrounds J/psi -> gamma eta_c (eta_c -> Xi- Xi+), " \
                 "J/psi -> pi- Lambda Sigma(1385)+ and " \
                 "J/psi -> Sigma(1385)- Sigma(1385)+ are found negligible after " \
                 "normalization to the total number of J/psi events.")

sel_jpsi_xi = build_tag_selection
sel_jpsi_xi.partial_rec([:Lambda, :pim]) {
              require_recoil_mass 1.312, 1.332    # anti-Xi+ signal region
            }

alg_jpsi_xi.apply(sel_jpsi_xi)
alg_jpsi_xi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_xi])

########################################################################
# Algorithm II: psi(3686) -> Xi- Xi+ (same final state and selection)
########################################################################
alg_psip_xi = Algorithm.new("PsipToXimXip")
alg_psip_xi.set_header(["PsipToXimXipAlg/PsipToXimXip.h"])
           .set_constant({ "ECMS" => [:double, 3.686] })
           .with_decay_card(decay_card_psip_xi)
           .note(:baryon_tag,
                 "Identical single baryon tag selection as J/psi -> Xi- Xi+: " \
                 "Xi- -> pi- Lambda, Lambda -> p pi-, antibaryon from the recoil mass " \
                 "of the pi-+ Lambda system.")
           .note(:psip_to_pipim_jpsi_veto,
                 "|M_recoil(pi+ pi-) - M_J/psi| > 0.005 GeV/c^2 to suppress " \
                 "psi(3686) -> pi+ pi- J/psi.")
           .note(:background_processes,
                 "Dominant backgrounds psi(3686) -> gamma chi_cJ with " \
                 "chi_cJ -> Xi- Xi+ and psi(3686) -> Sigma(1385)- Sigma(1385)+ are " \
                 "expected to populate the recoil-mass spectrum smoothly and are " \
                 "absorbed into the polynomial background shape of the fit.")
           .note(:angular_distribution,
                 "dN/d(cos theta) ~ 1 + alpha cos^2 theta fitted in 16 cos theta bins " \
                 "over [-0.8, 0.8].")

sel_psip_xi = build_tag_selection
sel_psip_xi.partial_rec([:Lambda, :pim]) {
              require_recoil_mass 1.308, 1.338    # anti-Xi+ signal region
            }

alg_psip_xi.apply(sel_psip_xi)
alg_psip_xi.execute_on([psip_data, psip_incMC, exMC_psip_xi])

########################################################################
# Algorithm III: J/psi -> Sigma(1385)-+ Sigma(1385)+- (single baryon tag)
########################################################################
alg_jpsi_sigma = Algorithm.new("JpsiToSigma1385mSigma1385p")
alg_jpsi_sigma.set_header(["JpsiToSigma1385mSigma1385pAlg/JpsiToSigma1385mSigma1385p.h"])
              .set_constant({ "ECMS" => [:double, 3.097] })
              .with_decay_card(decay_card_jpsi_sigma)
              .note(:baryon_tag,
                    "Single baryon tag: Sigma(1385)-+ -> pi-+ Lambda with " \
                    "Lambda -> p pi-; the anti-Sigma(1385)+- partner is inferred from " \
                    "the recoil mass of the pi-+ Lambda system. The charge-conjugate " \
                    "tag mode is covered by the conjugate multiplicity requirement.")
              .note(:tag_mass_window,
                    "The pi-+ Lambda system is required to satisfy " \
                    "|M(pi Lambda) - M_Sigma(1385)| < 0.035 GeV/c^2, with the " \
                    "candidate chosen by the minimum " \
                    "|M(pi Lambda) - M_Sigma(1385)(nominal)|. The charged " \
                    "Sigma(1385) mass width is 35-40 MeV (1 sigma level requirement); " \
                    "the window is varied by +-10 MeV/c^2 for the systematic uncertainty.")
              .note(:anti_baryon_signal_region,
                    "Signal events peak in the M_recoil(pi-+ Lambda) spectrum at the " \
                    "anti-Sigma(1385) nominal mass; the spectrum is fitted with the MC " \
                    "signal shape convolved with a Gaussian plus a second-order " \
                    "polynomial background. The peaking background in this mode is " \
                    "significant and is described by its own simulated shape with the " \
                    "normalized yield fixed in the fit.")
              .note(:angular_distribution,
                    "dN/d(cos theta) ~ 1 + alpha cos^2 theta fitted in 16 cos theta " \
                    "bins over [-0.8, 0.8].")
              .note(:background_processes,
                    "Backgrounds J/psi -> pi-+ Lambda anti-Sigma(1385)+-, " \
                    "J/psi -> Xi(1530)- anti-Xi+ + c.c. and " \
                    "J/psi -> Xi(1530)0 anti-Xi0 + c.c. are studied with generic MC.")

sel_jpsi_sigma = build_tag_selection
sel_jpsi_sigma.partial_rec([:Lambda, :pim]) {
                  require_recoil_mass 1.350, 1.420   # anti-Sigma(1385)+- region
                }

alg_jpsi_sigma.apply(sel_jpsi_sigma)
alg_jpsi_sigma.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_sigma])

########################################################################
# Algorithm IV: psi(3686) -> Sigma(1385)-+ Sigma(1385)+-
########################################################################
alg_psip_sigma = Algorithm.new("PsipToSigma1385mSigma1385p")
alg_psip_sigma.set_header(["PsipToSigma1385mSigma1385pAlg/PsipToSigma1385mSigma1385p.h"])
              .set_constant({ "ECMS" => [:double, 3.686] })
              .with_decay_card(decay_card_psip_sigma)
              .note(:baryon_tag,
                    "Identical single baryon tag selection as the " \
                    "J/psi -> Sigma(1385)-+ Sigma(1385)+- algorithm.")
              .note(:tag_mass_window,
                    "|M(pi Lambda) - M_Sigma(1385)| < 0.035 GeV/c^2; best candidate by " \
                    "the minimum |M(pi Lambda) - M_Sigma(1385)(nominal)|.")
              .note(:psip_to_pipim_jpsi_veto,
                    "The surviving backgrounds mainly come from " \
                    "psi(3686) -> pi+ pi- J/psi; the requirement " \
                    "|M_recoil(pi+ pi-) - M_J/psi| > 0.005 GeV/c^2 is applied.")
              .note(:angular_distribution,
                    "dN/d(cos theta) ~ 1 + alpha cos^2 theta fitted in 8 cos theta " \
                    "bins over [-0.8, 0.8].")

sel_psip_sigma = build_tag_selection
sel_psip_sigma.partial_rec([:Lambda, :pim]) {
                  require_recoil_mass 1.350, 1.420   # anti-Sigma(1385)+- region
                }

alg_psip_sigma.apply(sel_psip_sigma)
alg_psip_sigma.execute_on([psip_data, psip_incMC, exMC_psip_sigma])
