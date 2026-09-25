# ============================================================================
# Measurement of CP-even fractions F_+ for D0 -> pi+pi-pi0 and D0 -> K+K-pi0
# at psi(3770) using double-tag method with ~12 tag modes
# 2409.07197v1  —  BESIII
# ============================================================================

psip_data = DatasetManager.real_data.find("712_3773")
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================================
# Mode I: D0 -> pi+ pi- pi0
# ============================================================================

decay_card_pipipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_pi_pi_pi0"
  config.related_dataset = psip_data
  config.events = 200_000
  config.decay_card = decay_card_pipipi0
  config.cross_section = :default
end

alg_pipipi0 = Algorithm.new("D0ToPiPiPi0")
alg_pipipi0.set_header(["D0ToPiPiPi0Alg/D0ToPiPiPi0.h"])
            .set_constant({ "ECMS" => [:double, 3.773] })
            .set_alias({ "std::vector<double>" => "Vdouble" })

sel_pipipi0 = Selection.new
sel_pipipi0.select_track {
              cos_theta 0.93
              Vz 10.0
              Vxy 1.0
              nChrp ">=1"
              nChrn ">=1"
              nNet "==0"
            }
            .select_photon {
              tdc_emc_start 0
              tdc_emc_end 700
              angle_to_track 10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam ">=2"
            }
            .pid(method: :probability) {
              prob_cut 0.001
              identify :pion, against: [:kaon]
              npip ">=1"
              npim ">=1"
            }
            .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 25
              npi0 ">=1"
            }
            .kinematic_fit([:pip, :pim, :pi0]) {
              nominal
              constrain_four_momentum
              chi2_cut 200
            }

alg_pipipi0
  .note(:tag_modes, "12 tag modes: K+K-, pi+pi-, Kspipi, Kpi, Kspi0, Kspi0pi0, Ksomega, Kseta_gg, Ksetap_gammarho, Ksetap_pipieta, KL0pi0, KL0pi0pi0, KL0omega, K+K-pi0, pi+pi-pi0, pi+pi-pi+pi-, KSpipi (binned), KL0pipi (binned)")
  .note(:st_selection, "ST D0 candidates selected via MBC fit; DeltaE within +/-3 sigma per tag mode; |DeltaE| minimized for multiple candidates; MBC in [1.86, 1.87] for KL0 tags")
  .note(:dt_selection, "DT yield from 2D unbinned MBC(tag) vs MBC(sig) ML fit; signal described by 2D simulated shape convolved with Gaussian resolution")
  .note(:ks0_veto, "K_S0 veto on pi+pi-pi0 mode: reject if any pi+pi- combination has invariant mass in [0.481, 0.514] to suppress K_S0pi0 background")
  .note(:pid_details, "Kaon/pion PID using dE/dx and TOF likelihoods: L(K) > L(pi) for kaon, L(pi) > L(K) for pion")
  .note(:pi0_selection, "Photon pairs with gamma-gamma invariant mass in [0.115, 0.150]; both photons in barrel EMC region; kinematic fit constraining to pi0 mass")
  .note(:ks0_selection, "K_S0: two oppositely charged tracks |Vz|<20cm, vertex fit, mass [0.487,0.511], decay length > 2 sigma")
  .note(:eta_selection, "eta: gamma-gamma mass [0.505,0.575]; or pi+pi-pi0 mass [0.530,0.565]")
  .note(:omega_selection, "omega: pi+pi-pi0 mass [0.750,0.820]")
  .note(:rho0_selection, "rho0: pi+pi- mass [0.626,0.924]")
  .note(:etap_selection, "etap: pi+pi-eta(gg) mass [0.940,0.976]; or gamma rho0 mass [0.940,0.970]")
  .note(:cosmic_veto, "D->K+K- and D->pi+pi- modes: TOF flight time difference < 5 ns; Bhabha/dimuon PID veto")
  .note(:extra_track_veto, "D->K+K- and D->pi+pi- require >=1 extra shower >50 MeV or >=1 extra charged track")
  .note(:kl0_dt, "KL0X tag modes: partially reconstructed via Mmiss^2 fit; signal modeled by MC shape convolved with Gaussian; combinatorial background by 2nd-order Chebyshev")
  .note(:cp_fraction_method, "F_+ extracted from N+ and N- averaged over CP-eigen tag modes via least-squares fit; global CP-mixed and binned CP-mixed tag modes provide additional constraints")
  .with_decay_card(decay_card_pipipi0)
  .apply(sel_pipipi0)

alg_pipipi0.execute_on([psip_data, psip_incMC, exMC_pipipi0])

# ============================================================================
# Mode II: D0 -> K+ K- pi0
# ============================================================================

decay_card_kkpi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K+ K- pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_kkpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_K_K_pi0"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card_kkpi0
  config.cross_section = :default
end

alg_kkpi0 = Algorithm.new("D0ToKKPi0")
alg_kkpi0.set_header(["D0ToKKPi0Alg/D0ToKKPi0.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .set_alias({ "std::vector<double>" => "Vdouble" })

sel_kkpi0 = Selection.new
sel_kkpi0.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=1"
            nChrn ">=1"
            nNet "==0"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            nkp ">=1"
            nkm ">=1"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
          .kinematic_fit([:kp, :km, :pi0]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_kkpi0
  .note(:tag_modes, "12 tag modes: see Mode I notes for full list")
  .note(:st_selection, "ST D0 candidates selected via MBC fit; DeltaE within +/-3 sigma per tag mode; |DeltaE| minimized for multiple candidates; MBC in [1.86, 1.87] for KL0 tags")
  .note(:dt_selection, "DT yield from 2D unbinned MBC(tag) vs MBC(sig) ML fit; signal described by 2D simulated shape convolved with Gaussian resolution")
  .note(:pid_details, "Kaon/pion PID using dE/dx and TOF likelihoods: L(K) > L(pi) for kaon, L(pi) > L(K) for pion")
  .note(:pi0_selection, "Photon pairs with gamma-gamma invariant mass in [0.115, 0.150]; both photons in barrel EMC region; kinematic fit constraining to pi0 mass")
  .note(:ks0_selection, "K_S0: two oppositely charged tracks |Vz|<20cm, vertex fit, mass [0.487,0.511], decay length > 2 sigma")
  .note(:eta_selection, "eta: gamma-gamma mass [0.505,0.575]; or pi+pi-pi0 mass [0.530,0.565]")
  .note(:omega_selection, "omega: pi+pi-pi0 mass [0.750,0.820]")
  .note(:rho0_selection, "rho0: pi+pi- mass [0.626,0.924]")
  .note(:etap_selection, "etap: pi+pi-eta(gg) mass [0.940,0.976]; or gamma rho0 mass [0.940,0.970]")
  .note(:cosmic_veto, "D->K+K- and D->pi+pi- modes: TOF flight time difference < 5 ns; Bhabha/dimuon PID veto")
  .note(:extra_track_veto, "D->K+K- and D->pi+pi- require >=1 extra shower >50 MeV or >=1 extra charged track")
  .note(:kl0_dt, "KL0X tag modes: partially reconstructed via Mmiss^2 fit; signal modeled by MC shape convolved with Gaussian; combinatorial background by 2nd-order Chebyshev")
  .note(:cp_fraction_method, "F_+ extracted from N+ and N- averaged over CP-eigen tag modes via least-squares fit; global CP-mixed and binned CP-mixed tag modes provide additional constraints")
  .with_decay_card(decay_card_kkpi0)
  .apply(sel_kkpi0)

alg_kkpi0.execute_on([psip_data, psip_incMC, exMC_kkpi0])