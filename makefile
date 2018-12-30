COMPILER = gfortran
F77_OPTIONS = -g -o0 -C -ffixed-line-length-132 -funroll-loops
INCLUDE_DIRS  = -I/usr/local/include
LD_OPTIONS_NETCDF  = /usr/local/lib/libnetcdff.so

SPATH = ./IBIS
OBJPATH = ./build
DPATH = ./Debug

HEADERS = $(wildcard $(SPATH)/*.h)
SRC = $(wildcard $(SPATH)/*.f)
SRCNAME = $(notdir $(SRC))
OBJS = $(patsubst %.f, $(OBJPATH)/%.o, $(SRCNAME))
EXENAME = IBIS

$(EXENAME): $(OBJS) $(HEADERS)
	@mkdir -p $(DPATH)
	$(COMPILER) $(OBJS) $(F77_OPTIONS) $(INCLUDE_DIRS) $(LD_OPTIONS_NETCDF) -o $(DPATH)/$(EXENAME)

$(OBJS): $(OBJPATH)/%.o: $(SPATH)/%.f $(HEADERS)
	@mkdir -p $(OBJPATH)
	$(COMPILER) $(F77_OPTIONS) $(INCLUDE_DIRS) -c $< -o $@

all:
	make $(EXENAME)

clean:
	rm -rf $(OBJPATH)/*.o $(DPATH)/$(EXENAME) **/*.o

run:
	@ cd $(DPATH); ./$(EXENAME) --ini_infile_=./input/ibis.infile --diag_infile_=./input/diag.infile --cld_mon_=./input/cld.mon.nc --deltat_mon_=./input/deltat.nc --prec_mon_=./input/prec.mon.nc --rh_mon_=./input/rh.mon.nc --temp_mon_=./input/temp.mon.nc --trange_mon_=./input/trange.mon.nc --wetd_mon_=./input/wetd.mon.nc --wspd_mon_=./input/wspd.mon.nc --soita_sand_=./input/soita.sand.nc --soita_clay_=./input/soita.clay.nc --vegtype_=./input/vegtype.nc --surta_=./input/surta.nc --topo_=./input/topo.nc --deltat_=./input/deltat.nc --params_can_=./params/params.can --params_hyd_=./params/params.hyd --params_soi_=./params/params.soi --params_veg_=./params/params.veg --out_yearly_aet_=./output/yearly/aet.nc --out_yearly_biomass_=./output/yearly/biomass.nc --out_yearly_co2fluxes_=./output/yearly/co2fluxes.nc --out_yearly_csoi_=./output/yearly/csoi.nc --out_yearly_disturbf_=./output/yearly/disturbf.nc --out_yearly_exist_=./output/yearly/exist.nc --out_yearly_fcover_=./output/yearly/fcover.nc --out_yearly_npp_=./output/yearly/npp.nc --out_yearly_nsoi_=./output/yearly/nsoi.nc --out_yearly_plai_=./output/yearly/plai.nc --out_yearly_runoff_=./output/yearly/runoff.nc --out_yearly_sens_=./output/yearly/sens.nc --out_yearly_tsoi_=./output/yearly/tsoi.nc --out_yearly_vegtype0_=./output/yearly/vegtype0.nc --out_yearly_wsoi_=./output/yearly/wsoi.nc --out_yearly_zcanopy_=./output/yearly/zcanopy.nc --out_yearly_sapfrac_=./output/yearly/sapfrac.nc --out_yearly_dummyv_=./output/yearly/dummyv.nc --out_yearly_solar_=./output/yearly/solar.nc --out_yearly_albedo_=./output/yearly/albedo.nc --out_yearly_latent_=./output/yearly/latent.nc --out_yearly_totfall_=./output/yearly/totfall.nc --out_yearly_clitw_=./output/yearly/clitw.nc --out_yearly_csoislo_=./output/yearly/csoislo.nc --out_yearly_csoipas_=./output/yearly/csoipas.nc --out_monthly_aet_=./output/monthly/aet.nc --out_monthly_cloud_=./output/monthly/cloud.nc --out_monthly_co2ratio_=./output/monthly/co2ratio.nc --out_monthly_ir_=./output/monthly/ir.nc --out_monthly_lai_=./output/monthly/lai.nc --out_monthly_latent_=./output/monthly/latent.nc --out_monthly_npptot_=./output/monthly/npptot.nc --out_monthly_qa_=./output/monthly/qa.nc --out_monthly_rain_=./output/monthly/rain.nc --out_monthly_rh_=./output/monthly/rh.nc --out_monthly_runoff_=./output/monthly/runoff.nc --out_monthly_sens_=./output/monthly/sens.nc --out_monthly_snod_=./output/monthly/snod.nc --out_monthly_snof_=./output/monthly/snof.nc --out_monthly_snow_=./output/monthly/snow.nc --out_monthly_solar_=./output/monthly/solar.nc --out_monthly_temp_=./output/monthly/temp.nc --out_monthly_tsoi_=./output/monthly/tsoi.nc --out_monthly_wsoi_=./output/monthly/wsoi.nc --out_monthly_albedo_=./output/monthly/albedo.nc --out_monthly_dummyv_=./output/monthly/dummyv.nc --out_daily_rain_=./output/daily/rain.nc --out_daily_cloud_=./output/daily/cloud.nc --out_daily_rh_=./output/daily/rh.nc --out_daily_snow_=./output/daily/snow.nc --out_daily_aet_=./output/daily/aet.nc --out_daily_trunoff_=./output/daily/trunoff.nc --out_daily_srunoff_=./output/daily/srunoff.nc --out_daily_drainage_=./output/daily/drainage.nc --out_daily_wsoi_=./output/daily/wsoi.nc --out_daily_wisoi_=./output/daily/wisoi.nc --out_daily_snod_=./output/daily/snod.nc --out_daily_snof_=./output/daily/snof.nc --out_daily_co2ratio_=./output/daily/co2ratio.nc --out_daily_co2mic_=./output/daily/co2mic.nc --out_daily_templ_=./output/daily/templ.nc --out_daily_zcanopy_=./output/daily/zcanopy.nc --out_daily_laicanopy_=./output/daily/laicanopy.nc --out_global_=./output/ibis.out.global --out_vegtype_=./output/ibis.out.vegtype --out_yearsrun_=./output/ibis.out.yearsrun --out_diag_0_=./output/diag/0 --out_diag_1_=./output/diag/1 --out_diag_2_=./output/diag/2 --out_diag_3_=./output/diag/3 --out_diag_4_=./output/diag/4 --out_diag_5_=./output/diag/5 --out_diag_6_=./output/diag/6 --out_diag_7_=./output/diag/7 --out_diag_8_=./output/diag/8 --out_diag_9_=./output/diag/9 --temp_danom_= --trange_danom_= --prec_danom_= --cld_danom_= --rh_danom_= --wspd_danom_= --wetd_danom_= --prec_fanom_= --sphum_fanom_= --trange_fanomc_= --wspd_fanomc_=
