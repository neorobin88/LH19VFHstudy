rivet-mkhtml --mc-errs -s -c plots.conf -o Comparison/VBFvsGGH \
	     -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_pthj1 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthj1 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthj1_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_delta_r_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_pthj1 -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_pthj1_log -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_ht_log -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_pthjj12 -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_delta_y_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_m_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_delta_r_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_ht_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_delta_phi_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_delta_y_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_m_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_delta_r_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pth_log -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_ht -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_delta_r_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_pthj1_log -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_delta_phi_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_ht -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_delta_y_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_njets -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_delta_y_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_m_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthigh -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthlarge -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_delta_phi_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_m_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_njets -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pth -m /MC_HJETSVBF/atlas_pth -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_m_jj12 -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_delta_y_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_delta_r_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_ht -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthjj12 -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_m_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_pthjj12_log -m /MC_HJETSVBF/atlas_pth_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_ht_log -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_delta_y_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_pthjj12_log -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_pthjj12 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_njets -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_delta_r_jj12_log -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthjj12_log -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pthfull -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_pt2_pt1 -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_pt2_pt1 -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_pt2_pt1 -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_y_star -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_y_star -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_y_star -m /MC_HJETSVBF/rstudy_dr04_pth000_dy10_z_star -m /MC_HJETSVBF/rstudy_dr04_pth200_dy10_z_star -m /MC_HJETSVBF/rstudy_dr04_pth500_dy10_z_star \
	     --no-ratio \
	     yodafiles-v3/VBF-SH.yoda.gz:Title=VBF~MCNLO \
	     yodafiles-v3/VH-SH.yoda.gz:Title=ZH~MCNLO \
	     yodafiles-v3/HJJ-SH.yoda.gz:Title=HJJ~MCNLO \
	     yodafiles-v3/GGH-HEFT-SH.yoda.gz:Title=GGH~HEFT~MEPSNLO \
	     yodafiles-v3/GGH-SM-SH.yoda.gz:Title=GGH~SM~MEPSNLO
rivet-mkhtml --mc-errs -s -c plots.conf -o Comparison/VBFvsGGH2 -m .*MC_HJETSVBF.*dr[01][047].* -M .*MC_HJETSVBF/vbfvh.*dr[01][047].* \
	     --no-ratio \
	     yodafiles-v3/VBF-SH.yoda.gz:Title=VBF~MCNLO \
	     yodafiles-v3/VH-SH.yoda.gz:Title=ZH~MCNLO \
	     yodafiles-v3/HJJ-SH.yoda.gz:Title=HJJ~MCNLO \
	     yodafiles-v3/GGH-HEFT-SH.yoda.gz:Title=GGH~HEFT~MEPSNLO \
	     yodafiles-v3/GGH-SM-SH.yoda.gz:Title=GGH~SM~MEPSNLO