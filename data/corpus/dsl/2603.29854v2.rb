# e+e- -> K+ K- energy scan around psi(2S) (BESIII, 495 pb^-1, nine c.m. energies 3.58-3.71 GeV)
# Cross-section line-shape measurement with ConExc (mode 45 = K+ K-).

### Datasets ###
data_scan_1 = DatasetManager.real_data.find("704_psip_scan_1")   # 3581.5 MeV
data_scan_2 = DatasetManager.real_data.find("704_psip_scan_2")   # 3670.2 MeV
data_scan_3 = DatasetManager.real_data.find("704_psip_scan_3")   # 3680.1 MeV
data_scan_4 = DatasetManager.real_data.find("704_psip_scan_4")   # 3682.8 MeV
data_scan_5 = DatasetManager.real_data.find("704_psip_scan_5")   # 3684.2 MeV
data_scan_6 = DatasetManager.real_data.find("704_psip_scan_6")   # 3685.3 MeV
data_scan_7 = DatasetManager.real_data.find("704_psip_scan_7")   # 3686.5 MeV
data_scan_8 = DatasetManager.real_data.find("704_psip_scan_8")   # 3691.4 MeV
data_scan_9 = DatasetManager.real_data.find("704_psip_scan_9")   # 3709.8 MeV

data_points = [data_scan_1, data_scan_2, data_scan_3, data_scan_4, data_scan_5,
               data_scan_6, data_scan_7, data_scan_8, data_scan_9]

# ConExc decay card for e+e- -> K+ K- (mode 45); Particle vpho is auto-injected per energy point
decay_card_KK = <<~DECAYCARD
    Decay vpho
    1 ConExc 45;
    Enddecay

    Decay vhdr
    1 K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Signal MC over all nine scan energy points
exMC_signal_scan = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "psip_scan_KK_signal"
  config.events        = 100000
  config.decay_card    = decay_card_KK
  config.cross_section = :default
end
exMC_signal_scan.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "PsipScanKK"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.6865] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

event_selection = Selection.new
event_selection.select_track {
                  cos_theta  0.93
                  Vz         10.0
                  Vr         1.0
                  nChrp      "==1"
                  nChrn      "==1"
                  nNet       "==0"
                }
                .pid(method: :probability) {
                  prob_cut   0.001
                  identify :kaon, against: [:pion]
                  nkp        "==1"
                  nkm        "==1"
                }
                .kinematic_fit([:kp, :km]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

alg.note(:bhabha_veto,
         "Charged tracks must satisfy E/p < 0.8c to suppress Bhabha events (ratio of EMC deposited energy to MDC momentum).")
   .note(:cosmic_veto,
         "Cosmic ray veto: |T(K+) - T(K-)| < 3 ns using TOF flight times.")
   .note(:residual_momentum_cut,
         "Extra-track suppression: |p_CMS - p(K+) - p(K-)| < 0.07 GeV/c.")
   .note(:mkk_signal_region,
         "Signal region defined as M(K+K-) within 0.12 GeV around sqrt(s); applied in ROOT for unbinned maximum likelihood fit.")
   .note(:cross_section_analysis,
         "Cross-section line-shape fit is performed in the ROOT stage; ISR + beam-energy-spread convolution applied there.")

alg.with_decay_card(decay_card_KK).apply(event_selection)
alg.execute_on(data_points + exMC_signal_scan)
