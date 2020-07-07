function [EEE,EEE_HZ,ZZZ,calc_err,calc_err_edi,percent_err_replaced]=M3_err_flr_type_2(EEE,ZZZ,EEE_HZ,T,H,err_flr,err_flr2,err_flr3)
%Function for use with M3DET_v10.
%DC Edit June, 17, 2015
%Ensure function call has EEE_HZ as an output otherwise tipper error will
%not be rewritten

%DC November 18, 2014
%Calculates error percentage using Type 2 (see M3DET). Type 2 means that
%error percent is defined in two ways:
%For Zxx and Zxy components: Error_Percent = Error / abs(Zxy)
%For Zyy and Zyx components: Error_Percent = Error / abs(Zyx)
%
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
            if ic==1 || ic==2
                calc_err_edi(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(2,:,:));
            else
                calc_err_edi(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(3,:,:));
            end
        end
        
        N_nan = length(find(isnan(ZZZ(:))==1)); %Number of nan data points
        N = length(ZZZ(:))-N_nan; %Total number of non-nan data points
        count = 0; tip_count = 0;
        
        for ic=1:4
            for itrpo=1:length(H.d.x)
                for iu=1:nfim   %number of frequencies in *mat file
                    if ic==1
                        if  calc_err_edi(ic,iu,itrpo)<err_flr
                            count = count + 1;
                            EEE(ic,iu,itrpo)=err_flr*abs(ZZZ(2,iu,itrpo));   
                        end
                    elseif ic==2
                        if calc_err_edi(ic,iu,itrpo)<err_flr2
                            count = count + 1;
                            EEE(ic,iu,itrpo)=err_flr2*abs(ZZZ(2,iu,itrpo));
                        end
                    elseif ic==3
                        if  calc_err_edi(ic,iu,itrpo)<err_flr2
                            count = count + 1;
                            EEE(ic,iu,itrpo)=err_flr2*abs(ZZZ(3,iu,itrpo));   
                        end
                    elseif ic==4
                        if  calc_err_edi(ic,iu,itrpo)<err_flr
                            count = count + 1;
                            EEE(ic,iu,itrpo)=err_flr*abs(ZZZ(3,iu,itrpo));   
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
        %Perform backwards calculation to SOLVE FOR ERROR for errm3.
        %Identical conditional as above with ABS included. Simply solve for error for m3 data and edi
        %data
        for ic=1:4
            if ic==1 || ic==2
                calc_err(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(2,:,:));
            else
                calc_err(ic,:,:)=(EEE(ic,:,:))./abs(ZZZ(3,:,:));
            end
        end
        
        percent_err_replaced = count/N;
        percent_nan = N_nan/N;
        
end