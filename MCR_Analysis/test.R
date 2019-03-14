allquery <- amr.list.raw[unlist(lapply(amr.list.raw,
                                    function(x){sum(c("aac(3)-IV","blaCTX-M-14","blaEC","bleO","cmlA1","fosA3","qacL","sul2","sul3") %in% x) == 9}))]

mcrquery <- mcr.subset[unlist(lapply(mcr.subset,
                                    function(x){sum(c("aac(3)-IV","blaCTX-M-14","blaEC","bleO","cmlA1","fosA3","qacL","sul2","sul3", "mcr") %in% x) == 10}))]


matches <- mcr_rules_df %>%
  filter(all(lhs %in% unlist(strsplit("aph(3'')-Ib,aph(6)-Id,blaEC,bleO,clpK,oqxA,oqxB", split=","))) |
          all(lhs %in% unlist(strsplit("blaEC,floR,qnrS1", split=","))) |
          all(lhs %in% unlist(strsplit("aph(3')-Ia,blaEC,floR,tet(A)", split=","))) |
           all(lhs %in% unlist(strsplit("floR", split=","))) |
           all(lhs %in% unlist(strsplit("oqxB", split=","))) |
           all(lhs %in% unlist(strsplit("aadA1,aadA2,blaEC,cmlA1,dfrA12,qacL,sul3", split=","))) |
           all(lhs %in% unlist(strsplit("aac(3)-IIa,aadA1,aph(3'')-Ib,aph(6)-Id,blaEC,dfrA1,qacEdelta1,sul1,tet(A)", split=","))) |
           all(lhs %in% unlist(strsplit("bleO", split=","))) |
           all(lhs %in% unlist(strsplit("blaEC,sul3,tet(A)", split=","))) |
           all(lhs %in% unlist(strsplit("qnrS2", split=","))) |
           all(lhs %in% unlist(strsplit("sul3", split=","))) |
           all(lhs %in% unlist(strsplit("fosA3", split=","))) |
           all(lhs %in% unlist(strsplit("aac(3)-IId,aadA5,aph(3'')-Ib,aph(6)-Id,blaCTX-M-65,blaEC,blaTEM-1,dfrA17,mph(A),qacEdelta1,sul1,sul2,tet(A)", split=","))) |
           all(lhs %in% unlist(strsplit("oqxB2", split=","))))
                               
matches <- mcr_rules_df[(mcr_rules_df$lhs %in% unlist(strsplit("aph(3'')-Ib,aph(6)-Id,blaEC,bleO,clpK,oqxA,oqxB", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("blaEC,floR,qnrS1", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("aph(3')-Ia,blaEC,floR,tet(A)", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("floR", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("oqxB", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("aadA1,aadA2,blaEC,cmlA1,dfrA12,qacL,sul3", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("aac(3)-IIa,aadA1,aph(3'')-Ib,aph(6)-Id,blaEC,dfrA1,qacEdelta1,sul1,tet(A)", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("bleO", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("blaEC,sul3,tet(A)", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("qnrS2", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("sul3", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("fosA3", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("aac(3)-IId,aadA5,aph(3'')-Ib,aph(6)-Id,blaCTX-M-65,blaEC,blaTEM-1,dfrA17,mph(A),qacEdelta1,sul1,sul2,tet(A)", split=","))) |
                          (mcr_rules_df$lhs %in% unlist(strsplit("oqxB2", split=","))) | #####
                          all(unlist(strsplit("aph(3'')-Ib,aph(6)-Id,blaEC,bleO,clpK,oqxA,oqxB", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("blaEC,floR,qnrS1", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("aph(3')-Ia,blaEC,floR,tet(A)", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("floR", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("oqxB", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("aadA1,aadA2,blaEC,cmlA1,dfrA12,qacL,sul3", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("aac(3)-IIa,aadA1,aph(3'')-Ib,aph(6)-Id,blaEC,dfrA1,qacEdelta1,sul1,tet(A)", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("bleO", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("blaEC,sul3,tet(A)", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("qnrS2", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("sul3", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("fosA3", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("aac(3)-IId,aadA5,aph(3'')-Ib,aph(6)-Id,blaCTX-M-65,blaEC,blaTEM-1,dfrA17,mph(A),qacEdelta1,sul1,sul2,tet(A)", split=",")) %in% mcr_rules_df$lhs) |
                          all(unlist(strsplit("oqxB2", split=",")) %in% mcr_rules_df$lhs),]


matches <- mcr_rules_df[mcr_rules_df$lhs %in% unlist(strsplit("aph(3'')-Ib,aph(6)-Id,blaEC,bleO,clpK,oqxA,oqxB", split=",")),]
