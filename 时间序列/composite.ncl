begin
    fi = addfile("1960-2019-MAM-avemon-200HGT.nc","r")
    fe = addfile("1960-2019-MAM-avemon-U.nc","r")
    fd = addfile("1960-2019-MAM-avemon-V.nc","r")
    fa = addfile("1960-2019-MAM-avemon-PRE.nc","r")
    f1 = addfile("A1_projection.nc","r")
    f2 = addfile("A2_projection.nc","r")
    f3 = addfile("A3_projection.nc","r")
    A11 =f1->A1(:)
    A22 =f2->A2(:)
    A33 =f3->A3(:)
    A1_NOAA =dim_standardize(A11, 0)
    A2_NOAA =dim_standardize(A22, 0)
    A3_NOAA =dim_standardize(A33, 0)
    A1=abs(A1_NOAA)
    A2=abs(A2_NOAA)
    A3=abs(A3_NOAA)
    a=ind(A1.gt.A2.and.A1.gt.A3)
    A1_ind =A1_NOAA(a)
    A1_i=A1(a)
    b=ind(A1_ind.eq.A1_i)
    c=ind(A1_ind.ne.A1_i)

    z_z =fi->z(3:170,:,:)
    h_t =z_z/9.8
    h_t!0 ="time"
    h_t&time =z_z&time
    h_t!1 ="latitude"
    h_t&latitude =z_z&latitude
    h_t!2 ="longitude"
    h_t&longitude =z_z&longitude
    dim_h_t=dimsizes(h_t)
    h_t_mon=reshape(h_t,(/dim_h_t(0)/3,3,dim_h_t(1),dim_h_t(2)/))
    h1 = dim_avg_n_Wrap(h_t_mon, 1)   ;MAM平均
    h1!0 ="time"
    h1!1 ="latitude"
    h1&latitude =h_t&latitude
    h1!2 ="longitude"
    h1&longitude =h_t&longitude
    printVarSummary(h1)
    h_ano =dim_rmvmean_n_Wrap(h1, 0)
    printVarSummary(h_ano)
    h_ano_A1=h_ano(a,:,:)
    printVarSummary(h_ano_A1)
    h_ano_A1_pos=h_ano_A1(b,:,:)
    printVarSummary(h_ano_A1_pos)
    h_ano_A1_neg=h_ano_A1(c,:,:)
    printVarSummary(h_ano_A1_neg)
    h_ano_high =new((/18,721,1440/), float)
    h_ano_high(0:8,:,:)=h_ano_A1_pos
    h_ano_high(9:17,:,:)=-h_ano_A1_neg
    printVarSummary(h_ano_high)
    h_ano_ave =dim_avg_n_Wrap(h_ano_high, 0)
    printVarSummary(h_ano_ave)

    u = fe->u(3:170,{200},:,:)
    dim_u=dimsizes(u)
    u_mon=reshape(u,(/dim_u(0)/3,3,dim_u(1),dim_u(2)/))
    u1 =dim_avg_n_Wrap(u_mon, 1)  ;MAM平均
    u1!0 ="time"
    u1!1 ="latitude"
    u1&latitude =u&latitude
    u1!2 ="longitude"
    u1&longitude =u&longitude
    u_ano =dim_rmvmean_n_Wrap(u1, 0)
    printVarSummary(u_ano)
    u_ano_A1=u_ano(a,:,:)
    printVarSummary(u_ano_A1)
    u_ano_A1_pos=u_ano_A1(b,:,:)
    printVarSummary(u_ano_A1_pos)
    u_ano_A1_neg=u_ano_A1(c,:,:)
    printVarSummary(u_ano_A1_neg)
    u_ano_high =new((/18,721,1440/), float)
    u_ano_high(0:8,:,:)=u_ano_A1_pos
    u_ano_high(9:17,:,:)=-u_ano_A1_neg
    printVarSummary(u_ano_high)
    u_ano_ave =dim_avg_n_Wrap(u_ano_high, 0)
    printVarSummary(u_ano_ave)

    v = fd->v(3:170,{200},:,:)
    dim_v=dimsizes(v)
    v_mon=reshape(v,(/dim_v(0)/3,3,dim_v(1),dim_v(2)/))
    v1 =dim_avg_n_Wrap(v_mon, 1)  ;MAM平均
    v1!0 ="time"
    v1!1 ="latitude"
    v1&latitude =v&latitude
    v1!2 ="longitude"
    v1&longitude =v&longitude
    v_ano =dim_rmvmean_n_Wrap(v1, 0)
    printVarSummary(v_ano)
    v_ano_A1=v_ano(a,:,:)
    printVarSummary(v_ano_A1)
    v_ano_A1_pos=v_ano_A1(b,:,:)
    printVarSummary(v_ano_A1_pos)
    v_ano_A1_neg=v_ano_A1(c,:,:)
    printVarSummary(v_ano_A1_neg)
    v_ano_high =new((/18,721,1440/), float)
    v_ano_high(0:8,:,:)=v_ano_A1_pos
    v_ano_high(9:17,:,:)=-v_ano_A1_neg
    printVarSummary(v_ano_high)
    v_ano_ave =dim_avg_n_Wrap(v_ano_high, 0)
    printVarSummary(v_ano_ave)

    tp =fa->tp(3:170,:,:)
    pre = tp*1000
    pre!0 ="time"
    pre&time =tp&time
    pre!1 ="latitude"
    pre&latitude =tp&latitude
    pre!2 ="longitude"
    pre&longitude =tp&longitude
    dim_pre=dimsizes(pre)
    pre_mon=reshape(pre,(/dim_pre(0)/3,3,dim_pre(1),dim_pre(2)/))
    pre1 = dim_avg_n_Wrap(pre_mon, 1)  ;MAM平均
    pre1!0 ="time"
    pre1!1 ="latitude"
    pre1&latitude =pre&latitude
    pre1!2 ="longitude"
    pre1&longitude =v&longitude
    pre_ano =dim_rmvmean_n_Wrap(pre1, 0)
    printVarSummary(pre_ano)
    pre_ano_A1=pre_ano(a,:,:)
    printVarSummary(pre_ano_A1)
    pre_ano_A1_pos=pre_ano_A1(b,:,:)
    printVarSummary(pre_ano_A1_pos)
    pre_ano_A1_neg=pre_ano_A1(c,:,:)
    printVarSummary(pre_ano_A1_neg)
    pre_ano_high =new((/18,721,1440/), float)
    pre_ano_high(0:8,:,:)=pre_ano_A1_pos
    pre_ano_high(9:17,:,:)=-pre_ano_A1_neg
    printVarSummary(pre_ano_high)
    pre_ano_ave =dim_avg_n_Wrap(pre_ano_high, 0)
    printVarSummary(pre_ano_ave)

     ;高指数的显著性检验
     n_high       =dimsizes(pre_ano_high(:,0,0))               
     df_high      = n_high-1   
     ;ts显著性检验
     pre_comp   =dim_avg_n(pre_ano_high(:,:,:),0)
     pre_std=dim_variance_n_Wrap(pre_ano_high(:,:,:),0)
     pre_std=sqrt(pre_std/n_high)
     pre_std=where(pre_std.eq.0,pre_std@_FillValue,pre_std)
     t_pre=pre_comp/pre_std
     prob_ts_high=student_t(t_pre, df_high)
     copy_VarMeta(pre_ano_ave,prob_ts_high)
     printVarSummary(prob_ts_high)

    f_out_name = "com_A1_200_world.nc"
    system("rm -f " + f_out_name)
    f_out = addfile(f_out_name, "c")
    f_out->u_ano_ave = u_ano_ave
    f_out->v_ano_ave = v_ano_ave
    f_out->h_ano_ave = h_ano_ave
    f_out->pre_ano_ave = pre_ano_ave
    f_out->prob_ts_high = prob_ts_high
end