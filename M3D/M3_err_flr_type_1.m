function [EEE,EEE_HZ,ZZZ,calc_err,calc_err_edi,percent_err_replaced]=M3_err_flr_type_1(EEE,ZZZ,EEE_HZ,T,H,err_flr,err_flr2,err_flr3)
%Function for use with M3DET_v10.
%DC Edit June, 17, 2015
%Ensure function call has EEE_HZ as an output otherwise tipper error will
%not be rewritten

%DC November 18, 2014
%Calculates error percentage using Type 1 (see M3DET). Type 1 means that
%error percent is defined in one way for all Z components:
%Error_Percent = Error / sqrt(abs(Zxy*Zyx))

%Syntax:
%[EEE,calc_err,calc_err_edi]=err_flr_type_2(EEE,ZZZ,EEE_HZ,T,H,err_flr,err_flr2)
%Inputs:
%EEE=EDI error
%ZZZ=EDI impedances
%EEE_HZ=Tipper error
%T=period
%err_flr=error floor to be applied to Zxx and Zyy
%err_flr2=error floor to be applied to Zxy and Zyx
%H=gui H

    nfim=length(T);
    calc_err_edi=[];calc_err=[];
    for ic=1:4
        calc_err_edi(ic,:,:)=(real(EEE(ic,:,:)))./sqrt(abs(ZZZ(2,:,:).*(ZZZ(3,:,:))));
    end
    
    N = length(ZZZ(:))-length(find(isnan(ZZZ(:))==1)); %Total number of non-nan impedance data points
    N_tip = length(EEE_HZ(:))-length(find(isnan(EEE_HZ(:))==1)); %Number of non-nan tipper data points
    xxcount = 0; xycount = 0; tip_count = 0;
    for ic=1:4
        for itrpo=1:length(H.d.x)
            for iu=1:nfim   %number of frequencies in *mat file
                if ic==1 || ic==4
                    if  calc_err_edi(ic,iu,itrpo)<err_flr
                        xycount = xycount + 1;
                        EEE(ic,iu,itrpo)=err_flr*sqrt(abs(ZZZ(2,iu,itrpo).*ZZZ(3,iu,itrpo)));   
                    end
                elseif ic==2 || ic==3
                    if  calc_err_edi(ic,iu,itrpo)<err_flr2
                        xxcount = xxcount + 1;
                        EEE(ic,iu,itrpo)=err_flr2*sqrt(abs(ZZZ(2,iu,itrpo).*ZZZ(3,iu,itrpo)));   
                    end
                end
                if (ic ==1 || ic == 2) && (H.resp == 2 || H.resp == 4) %tipper
                    if EEE_HZ(ic,iu,itrpo) < err_flr3
                        tip_count = tip_count + 1;
                        EEE_HZ(ic,iu,itrpo)=err_flr3;   %absolute error floor for tipper error
                    end
                end
            end
        end
    end 
    %%
    %Perform backwards calculation to SOLVE FOR ERROR for errm3.
    %Identical conditional as above with ABS included. Simply solve for error for m3 data and edi
    %data
    for ic=1:4
    calc_err(ic,:,:)=(real(EEE(ic,:,:)))./sqrt(abs(ZZZ(2,:,:).*(ZZZ(3,:,:))));
    end
    
    percent_err_replaced = (xxcount+xycount)/(N);
    
    
        
end