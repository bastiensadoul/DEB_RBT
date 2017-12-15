function Hleg = shrbt(nr)
  % makes graphs of survivor functions of various quantities
  % nr: vector with selected graphs
  % Hleg: hangle of legend
  
  close all

  llegend_rbt = {...
    {'-', 2, [1 0 0]}, 'Oncorhynchus_mykiss'; ....
    {'-', 2, [0 0 1]}, 'Salmoniformes'; ....
    {'-', 2, [0 0 0]}, 'Actinopterygii'; ....
  };
  shstat_options('default');
  shstat_options('y_label', 'on'); % if 'off' (default), no `survivor function' shown on yaxis

  for i = nr
  switch i
      
  case 1 % T_typical
  fig(i);
  shstat_options('x_transform', 'none');
  T_typical = K2C(read_allStat('T_typical'));
  [Hfig, Hleg, val] = shstat(T_typical, llegend_rbt,[], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_Bor = val(select_01('Boreogadus'));
  fprintf(['T_typ Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('typical temperature, {^\circ}C')

  case 2 % kappa
  fig(i);
  shstat_options('x_transform', 'none');
  [Hfig, Hleg, val] = shstat({'kap'}, llegend_rbt,[], i);
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['kap Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('allocation fraction \kappa, -')
 
  case 3 % [p_M]
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'p_M'}, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['[p_M] Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('spec somatic maintenance _{10}log [p_M], J/d.cm^3')

  case 4 % v
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('v', 's_M'); v = prod(vars,2);
  [Hfig, Hleg, val] = shstat(v, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['v Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('post-metam energy conductance _{10}log v, cm/d')

  case 5 % [E_m]
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'E_m'}, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['[E_m] Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('reserve capacity _{10}log [E_m], J/cm^3')

  case 6 % E_0
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'E_0'}, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['E_0 Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('initial energy _{10}log E_0, J')

  case 7 % Ww_i
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'Ww_i'}, llegend_rbt,[], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['Ww_i Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('ultimate body weight _{log} W_w^\infty, g')
  
  case 8 % r: dWm/W_dWm
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('dWm', 'W_dWm', 'c_T'); r = vars(:,1)./vars(:,2)./vars(:,3);
  [Hfig, Hleg, val] = shstat(r, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['r Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('max spec growth _{10}log r, 1/d')

  case 9 % R_i/Ww_i
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('R_i', 'c_T', 'Ww_i'); R_i = vars(:,1)./vars(:,2)./ vars(:,3);
  [Hfig, Hleg, val] = shstat(R_i, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['R_i/Ww_i Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('max reprod rate _{10}log R_\infty, #/d')
  
  case 10 % p_Ri/Ww_i
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('p_Ri', 'c_T', 'Ww_i'); p_Ri = vars(:,1)./vars(:,2)./vars(:,3);
  [Hfig, Hleg, val] = shstat(p_Ri, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['p_Ri/Ww_i Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('reproductive investment _{10}log p_R^\infty, J/d')

  case 11 % g
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'g'}, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['g Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('energy investment ratio, _{10}log g')

  case 12 % s_M
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig Hleg, val] = shstat({'s_M'}, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['s_M Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('accelration factor _{10}log s_M')

  case 13 % s_s
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'s_s'}, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['s_s Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('supply stress _{10}log s_s')

  case 14 % s_H^bp
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'s_Hbp'}, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['s_Hbp Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('maturity ratio _{10}log s_H^{bp}')

  case 15 % s_HL^bp
  fig(i);
  shstat_options('x_transform', 'none');
  [Hfig, Hleg, val] = shstat({'s_HLbp'}, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['s_HLbp Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('maturity density ratio s_{HL}^{bp}')

  case 16 % h_a
  fig(i);
  shstat_options('x_transform', 'log10');
  [Hfig, Hleg, val] = shstat({'h_a'}, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['h_a Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('aging acceleration _{10}log h_a, 1/d')

  case 17 % a_m
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('a_m', 'c_T'); a_m = vars(:,1).*vars(:,2);
  [Hfig, Hleg, val] = shstat(a_m, llegend_rbt, 'At T_{ref}', i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['a_m Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('longevity _{10}log a_m, d')

  case 18 % p_Mi/p_Ai
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('p_Si', 'p_Ai'); s_SA = vars(:,1)./vars(:,2);
  [Hfig, Hleg, val] = shstat(s_SA, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['p_Si/p_Ai Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('_{10}log p_M^\infty/ p_A^\infty')

  case 19 % p_Ji/p_Ai
  fig(i);
  shstat_options('x_transform', 'log10');
  vars = read_allStat('p_Ji', 'p_Ai'); s_JA = vars(:,1)./vars(:,2);
  [Hfig, Hleg, val] = shstat(s_JA, llegend_rbt, [], i);      
  med_Act = median(val(select_01('Actinopterygii'))); med_Sal = median(val(select_01('Gadiformes'))); med_rbt = median(val(select_01('Boreogadus')));
  fprintf(['p_Ji/p_Ai Act ', num2str(med_Act), '; Sal ', num2str(med_Sal), '; rbr ', num2str(med_rbt), '\n'])
  close(Hleg);
  fig(Hfig);
  xlabel('_{10}log p_J^\infty/ p_A^\infty')

  end
  end

  Hleg = shllegend(llegend_rbt);