function export_matlab_reference(source_dir, output_file)
% Run unmodified MATLAB modules from the pinned checkout in MATLAB or Octave.
addpath(source_dir);
Forc = struct('Temp_C',20,'Pres_kPa',101.325,'Wind_ms',2,'VapP_kPa',1.5, ...
    'SRad_Wm2',200,'LRad_Wm2',350,'Albd',0.2,'LAI',2,'P',1e-6,'fv',0.6);
Veg = struct('zr',10,'h_c',20,'d',13,'z0m',2,'z0h',0.2,'rc_min',100, ...
    'cw_per_LAI',0.0012,'free_throughfall_f0',0);
PM = hsbsm_pm(Forc,Veg,struct('rs_soil',157.70),struct());
pm = [PM.e_wetCan,PM.e_dryCan,PM.e_bare,PM.ra,PM.Rn];
[pn,ec,can] = hsbsm_canopy_pm(Forc,struct('S_can_m',0.0005),Veg,3600, ...
    struct('external_e_wetCan_mps',PM.e_wetCan));
canopy = [pn,ec,can.S_can_m];
zone = struct('z',0.65,'theta_min',0.045,'theta_sat',0.38,'theta0',0.30, ...
    'vG',struct('alpha',6,'n',2.2,'Ks',2e-6));
Soil = struct('RZ',zone,'TZ',zone);
Soil.TZ.z = 0.35; Soil.TZ.cap_frac_of_K = 0.3;
P = struct('Soil',Soil,'Infil',struct('model','Milly','K_mps',2e-6, ...
    'F0_m',0.001,'Fmin_m',0.01,'resetDrySteps',3,'max_c',0.7));
[inf,rex,~,si] = hsbsm_infil_milly([1e-5;1e-7;0], ...
    struct('Fi_m',[0;0.02;0.1],'dryCount',[0;0;3]),P,3600,[0.2;0.3;0.38]);
infiltration = [inf,rex,si.Fi_m,si.dryCount,si.last_InfCap];
P.Infil.model = 'Troch';
P.Troch.alpha_inv_m = 1/3;
[inf,rex,~,si] = hsbsm_infil_milly([1e-5;1e-7;0], ...
    struct('Fi_m',[0;0.02;0.1],'dryCount',[0;0;3]),P,3600,[0.2;0.3;0.38]);
troch_override = [inf,rex,si.Fi_m,si.dryCount,si.last_InfCap];
% Zero aquifer storage keeps the source unit bug inactive in this parity case.
[sr,st,fl] = hsbsm_rz_tz_coupled([1e-7;2e-7],1e-8,2e-9, ...
    struct('theta',[0.30;0.32]),struct('theta',[0.29;0.31]),Soil,60, ...
    struct('SAq',struct('h_m',[0;0],'S',[0;0]),'A_node_m2',[100;100]));
soil = [sr.theta,st.theta,fl.ET_rz_mps,fl.E_soil_mps,fl.q_rz2tz_mps,fl.q_tz2aq_mps];
x = [5;17;35;60]; w = [8;7;5;3];
ap = struct('k',1e-4,'f',0.2,'iota',0.1,'D',1,'S0',0.2*w.*[0.1;0.2;0.15;0.1]);
opts = struct('theta',1,'omega',0.8,'picard_max',60,'picard_tol',1e-11, ...
    'verbose',false,'mass_check',true,'plot_live',false,'safeguard',true);
[aq_q,aq_over,~,aq_h,aq_s,aq_diag] = hsB_solver(x,w,ones(10,1)*1e-7,60,10,ap,opts,false);
save(output_file,'pm','canopy','infiltration','troch_override','soil','aq_q','aq_over','aq_h','aq_s','-v7');
fprintf('Wrote reference fixture: %s\n',output_file);
end
