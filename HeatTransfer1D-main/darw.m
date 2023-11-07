
clearvars -except res1


plot(obsGST(NTB:NTB+4382,1),'k')
hold on
plot(res_ini(6,:),'b')

plot(res_68(6,:),'r')


plot(res_69(6,:),'g')



