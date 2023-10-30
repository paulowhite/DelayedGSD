kMax=2
k=2
z=-3
Info.d=12.58798
Info.i=c(10.60272, 20.74867)
ck =1.49773
ck.unrestricted =1.49773
lk=c(0.2433581, 1.9961649)
uk=c(2.545884, 1.996165)
reason.interim =c("no boundary crossed", NA)
method=1
bindingFutility=TRUE
cNotBelowFixedc=FALSE
continuity.correction=0

FinalPvalue(Info.d = Info.d[1:min(k,kMax-1)],
            Info.i = Info.i[1:k],
            ck = ck[1:min(k,kMax-1)], ck.unrestricted = ck.unrestricted[1:min(k,kMax-1)],
            lk = lk[1:k],
            uk = uk[1:k],
            kMax = kMax,
            delta = 0, 
            estimate = ifelse(k<kMax, z / sqrt(Info.d[k]), z / sqrt(Info.i[k])),
            reason.interim = reason.interim[1:k],
            method = method,
            bindingFutility = bindingFutility, 
            cNotBelowFixedc = cNotBelowFixedc,
            continuity.correction = continuity.correction)

Info.d = c(12.58797814)
Info.i = c(10.60271846, 20.74866926)
ck = c(1.959964,1.9961649)
ck.unrestricted = c(1.49772992,1.9961649)
lk = c(0.24335814, 1.9961649)
uk = c(2.5458844, 1.9961649)
kMax = 2
reason.interim = c("no boundary crossed",NA)
method = 1
bindingFutility = TRUE
cNotBelowFixedc = TRUE
continuity.correction = 0
xlim = c(0.7,2.2) 
title = "no correction"

k=2
z=-3

FinalPvalue(Info.d = Info.d[1:min(k,kMax-1)],
            Info.i = Info.i[1:k],
            ck = ck[1:min(k,kMax)], ck.unrestricted = ck.unrestricted[1:min(k,kMax)],
            lk = lk[1:k],
            uk = uk[1:k],
            kMax = kMax,
            delta = 0, 
            estimate = ifelse(k<kMax, z / sqrt(Info.d[k]), z / sqrt(Info.i[k])),
            reason.interim = reason.interim[1:k],
            method = method,
            bindingFutility = bindingFutility, 
            cNotBelowFixedc = cNotBelowFixedc,
            continuity.correction = continuity.correction)




Info.d = c(12.58797814)
Info.i = c(10.60271846, 20.74866926)
ck = c(1.959964,1.9961649)
ck.unrestricted = c(1.49772992,1.9961649)
lk = c(0.24335814, 1.9961649)
uk = c(2.5458844, 1.9961649)
kMax = 2
reason.interim = c("no boundary crossed",NA)
method = 1
bindingFutility = TRUE
cNotBelowFixedc = TRUE
continuity.correction = 1
xlim = c(0.7,2.3) 
title = "no correction"

k=1
z=1.95997

FinalPvalue(Info.d = Info.d[1:min(k,kMax-1)],
            Info.i = Info.i[1:k],
            ck = ck[1:min(k,kMax)], ck.unrestricted = ck.unrestricted[1:min(k,kMax)],
            lk = lk[1:k],
            uk = uk[1:k],
            kMax = kMax,
            delta = 0, 
            estimate = ifelse(k<kMax, z / sqrt(Info.d[k]), z / sqrt(Info.i[k])),
            reason.interim = reason.interim[1:k],
            method = method,
            bindingFutility = bindingFutility, 
            cNotBelowFixedc = cNotBelowFixedc,
            continuity.correction = continuity.correction)
