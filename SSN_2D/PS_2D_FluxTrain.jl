using Flux, Flux.Tracker
include("PS_Loss.jl")

#define model parameters
# Jee = param(0.5)
# Jei = param(-0.5)
# Jie = param(2.7)
# Jii = param(-1.5)

J0 = param([0.5 -0.5; 2.7 -1.5])
i2e = param(0.8)

#define idealSpect
idealSp = [0.0040    0.1023    0.0697    0.0483
    0.0040    0.0907    0.0625    0.0437
    0.0040    0.0764    0.0535    0.0378
    0.0040    0.0688    0.0486    0.0346
    0.0040    0.0649    0.0460    0.0329
    0.0040    0.0627    0.0445    0.0319
    0.0039    0.0614    0.0436    0.0312
    0.0039    0.0606    0.0429    0.0307
    0.0039    0.0600    0.0424    0.0303
    0.0039    0.0595    0.0420    0.0300
    0.0038    0.0592    0.0416    0.0296
    0.0038    0.0589    0.0412    0.0293
    0.0038    0.0586    0.0409    0.0290
    0.0037    0.0584    0.0406    0.0287
    0.0037    0.0583    0.0402    0.0284
    0.0037    0.0582    0.0399    0.0281
    0.0036    0.0581    0.0396    0.0278
    0.0036    0.0580    0.0393    0.0275
    0.0035    0.0580    0.0390    0.0272
    0.0035    0.0580    0.0387    0.0269
    0.0034    0.0580    0.0384    0.0265
    0.0034    0.0581    0.0382    0.0262
    0.0033    0.0582    0.0379    0.0259
    0.0033    0.0584    0.0377    0.0256
    0.0032    0.0586    0.0374    0.0253
    0.0032    0.0589    0.0372    0.0250
    0.0031    0.0592    0.0370    0.0247
    0.0031    0.0596    0.0368    0.0244
    0.0030    0.0601    0.0366    0.0242
    0.0030    0.0606    0.0364    0.0239
    0.0029    0.0612    0.0363    0.0236
    0.0029    0.0618    0.0362    0.0234
    0.0028    0.0626    0.0361    0.0231
    0.0028    0.0635    0.0360    0.0229
    0.0027    0.0644    0.0359    0.0226
    0.0027    0.0655    0.0359    0.0224
    0.0026    0.0667    0.0359    0.0222
    0.0026    0.0680    0.0359    0.0220
    0.0025    0.0694    0.0359    0.0218
    0.0024    0.0710    0.0360    0.0216
    0.0024    0.0728    0.0361    0.0214
    0.0023    0.0748    0.0362    0.0213
    0.0023    0.0769    0.0364    0.0211
    0.0022    0.0793    0.0366    0.0209
    0.0022    0.0820    0.0368    0.0208
    0.0022    0.0849    0.0371    0.0207
    0.0021    0.0881    0.0374    0.0206
    0.0021    0.0917    0.0378    0.0204
    0.0020    0.0956    0.0382    0.0204
    0.0020    0.0999    0.0386    0.0203
    0.0019    0.1048    0.0391    0.0202
    0.0019    0.1101    0.0397    0.0201
    0.0018    0.1160    0.0403    0.0201
    0.0018    0.1225    0.0410    0.0200
    0.0018    0.1297    0.0417    0.0200
    0.0017    0.1376    0.0426    0.0200
    0.0017    0.1463    0.0435    0.0200
    0.0016    0.1559    0.0445    0.0200
    0.0016    0.1662    0.0456    0.0200
    0.0016    0.1772    0.0468    0.0200
    0.0015    0.1889    0.0481    0.0201
    0.0015    0.2010    0.0496    0.0201
    0.0015    0.2130    0.0511    0.0202
    0.0014    0.2245    0.0528    0.0203
    0.0014    0.2348    0.0547    0.0204
    0.0014    0.2431    0.0568    0.0205
    0.0013    0.2485    0.0591    0.0206
    0.0013    0.2501    0.0615    0.0208
    0.0013    0.2476    0.0642    0.0209
    0.0013    0.2408    0.0672    0.0211
    0.0012    0.2301    0.0704    0.0213
    0.0012    0.2162    0.0740    0.0215
    0.0012    0.2000    0.0778    0.0217
    0.0011    0.1827    0.0820    0.0219
    0.0011    0.1652    0.0865    0.0222
    0.0011    0.1482    0.0914    0.0225
    0.0011    0.1322    0.0966    0.0227
    0.0010    0.1175    0.1021    0.0230
    0.0010    0.1043    0.1078    0.0233
    0.0010    0.0924    0.1137    0.0236
    0.0010    0.0819    0.1195    0.0239
    0.0010    0.0727    0.1252    0.0243
    0.0009    0.0647    0.1304    0.0246
    0.0009    0.0576    0.1348    0.0249
    0.0009    0.0514    0.1382    0.0253
    0.0009    0.0460    0.1401    0.0256
    0.0009    0.0413    0.1404    0.0259
    0.0008    0.0371    0.1388    0.0262
    0.0008    0.0335    0.1353    0.0265
    0.0008    0.0303    0.1302    0.0268
    0.0008    0.0274    0.1236    0.0270
    0.0008    0.0249    0.1159    0.0272
    0.0007    0.0227    0.1075    0.0274
    0.0007    0.0207    0.0988    0.0275
    0.0007    0.0190    0.0902    0.0275
    0.0007    0.0174    0.0818    0.0275
    0.0007    0.0160    0.0739    0.0274
    0.0007    0.0147    0.0666    0.0273
    0.0007    0.0135    0.0599    0.0270
    0.0006    0.0125    0.0539    0.0267
    0.0006    0.0116    0.0484    0.0263]

idealSp = idealSp./mean(idealSp)

ll, conv, spect = PS_Loss(J0, i2e, idealSp)

theta = Params([Jee, Jei, Jie, Jii, i2e])
grads = Tracker.gradient(() -> loss(Jee, Jei, Jie, Jii, i2e, idealSp), theta)

using Flux.Tracker: grad, update!

#I want to optimize using gradient descent, with a learning rate of 0.1
opt = Descent(0.1)

for p in (Jee, Jei, Jie, Jii, i2e)
    update!(opt, p, grads[p])
end
