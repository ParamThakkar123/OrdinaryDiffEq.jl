struct Rosenbrock23Tableau{T}
    c₃₂::T
    d::T
end

function Rosenbrock23Tableau(T)
    c₃₂ = convert(T, 6 + sqrt(2))
    d = convert(T, 1 / (2 + sqrt(2)))
    Rosenbrock23Tableau(c₃₂, d)
end

struct Rosenbrock32Tableau{T}
    c₃₂::T
    d::T
end

function Rosenbrock32Tableau(T)
    c₃₂ = convert(T, 6 + sqrt(2))
    d = convert(T, 1 / (2 + sqrt(2)))
    Rosenbrock32Tableau(c₃₂, d)
end

struct ROS3PTableau{T, T2}
    a21::T
    a31::T
    a32::T
    C21::T
    C31::T
    C32::T
    b1::T
    b2::T
    b3::T
    btilde1::T
    btilde2::T
    btilde3::T
    gamma::T2
    c2::T2
    c3::T2
    d1::T
    d2::T
    d3::T
end

function ROS3PTableau(T, T2)
    gamma = convert(T, 1 / 2 + sqrt(3) / 6)
    igamma = inv(gamma)
    a21 = convert(T, igamma)
    a31 = convert(T, igamma)
    a32 = convert(T, 0)
    C21 = convert(T, -igamma^2)
    tmp = -igamma * (convert(T, 2) - convert(T, 1 / 2) * igamma)
    C31 = -igamma * (convert(T, 1) - tmp)
    C32 = tmp
    tmp = igamma * (convert(T, 2 / 3) - convert(T, 1 / 6) * igamma)
    b1 = igamma * (convert(T, 1) + tmp)
    b2 = tmp
    b3 = convert(T, 1 / 3) * igamma
    # btilde1 = convert(T,2.113248654051871)
    # btilde2 = convert(T,1.000000000000000)
    # btilde3 = convert(T,0.4226497308103742)
    btilde1 = b1 - convert(T, 2.113248654051871)
    btilde2 = b2 - convert(T, 1.000000000000000)
    btilde3 = b3 - convert(T, 0.4226497308103742)
    c2 = convert(T, 1)
    c3 = convert(T, 1)
    d1 = convert(T, 0.7886751345948129)
    d2 = convert(T, -0.2113248654051871)
    d3 = convert(T, -1.077350269189626)
    ROS3PTableau(
        a21, a31, a32, C21, C31, C32, b1, b2, b3, btilde1, btilde2, btilde3, gamma,
        c2, c3, d1, d2, d3)
end

struct Rodas3Tableau{T, T2}
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    C21::T
    C31::T
    C32::T
    C41::T
    C42::T
    C43::T
    b1::T
    b2::T
    b3::T
    b4::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    gamma::T2
    c2::T2
    c3::T2
    d1::T
    d2::T
    d3::T
    d4::T
end

function Rodas3Tableau(T, T2)
    gamma = convert(T, 1 // 2)
    a21 = convert(T, 0)
    a31 = convert(T, 2)
    a32 = convert(T, 0)
    a41 = convert(T, 2)
    a42 = convert(T, 0)
    a43 = convert(T, 1)
    C21 = convert(T, 4)
    C31 = convert(T, 1)
    C32 = convert(T, -1)
    C41 = convert(T, 1)
    C42 = convert(T, -1)
    C43 = convert(T, -8 // 3)
    b1 = convert(T, 2)
    b2 = convert(T, 0)
    b3 = convert(T, 1)
    b4 = convert(T, 1)
    btilde1 = convert(T, 0.0)
    btilde2 = convert(T, 0.0)
    btilde3 = convert(T, 0.0)
    btilde4 = convert(T, 1.0)
    c2 = convert(T, 0.0)
    c3 = convert(T, 1.0)
    c4 = convert(T, 1.0)
    d1 = convert(T, 1 // 2)
    d2 = convert(T, 3 // 2)
    d3 = convert(T, 0)
    d4 = convert(T, 0)
    Rodas3Tableau(a21, a31, a32, a41, a42, a43, C21, C31, C32, C41, C42, C43, b1, b2, b3,
        b4, btilde1, btilde2, btilde3, btilde4, gamma, c2, c3, d1, d2, d3, d4)
end

struct Rodas3PTableau{T, T2}
    a21::T
    a41::T
    a42::T
    a43::T
    C21::T
    C31::T
    C32::T
    C41::T
    C42::T
    C43::T
    C51::T
    C52::T
    C53::T
    C54::T
    gamma::T
    c2::T2
    c3::T2
    d1::T
    d2::T
    d3::T
    h21::T
    h22::T
    h23::T
    h24::T
    h25::T
    h31::T
    h32::T
    h33::T
    h34::T
    h35::T
    h2_21::T
    h2_22::T
    h2_23::T
    h2_24::T
    h2_25::T
end

function Rodas3PTableau(T, T2)
    gamma = convert(T, 1 // 3)
    a21 = convert(T, 4.0 / 3.0)
    a41 = convert(T, 2.90625)
    a42 = convert(T, 3.375)
    a43 = convert(T, 0.40625)
    C21 = -convert(T, 4.0)
    C31 = convert(T, 8.25)
    C32 = convert(T, 6.75)
    C41 = convert(T, 1.21875)
    C42 = -convert(T, 5.0625)
    C43 = -convert(T, 1.96875)
    C51 = convert(T, 4.03125)
    C52 = -convert(T, 15.1875)
    C53 = -convert(T, 4.03125)
    C54 = convert(T, 6.0)
    c2 = convert(T2, 4.0 / 9.0)
    c3 = convert(T2, 0.0)
    d1 = convert(T, 1.0 / 3.0)
    d2 = -convert(T, 1.0 / 9.0)
    d3 = convert(T, 1.0)
    h21 = convert(T, 1.78125)
    h22 = convert(T, 6.75)
    h23 = convert(T, 0.15625)
    h24 = -convert(T, 6.0)
    h25 = -convert(T, 1.0)
    h31 = convert(T, 4.21875)
    h32 = -convert(T, 15.1875)
    h33 = -convert(T, 3.09375)
    h34 = convert(T, 9.0)
    h35 = convert(T, 0.0)
    h2_21 = convert(T, 4.21875)
    h2_22 = -convert(T, 2.025)
    h2_23 = -convert(T, 1.63125)
    h2_24 = -convert(T, 1.7)
    h2_25 = -convert(T, 0.1)
    Rodas3PTableau(a21, a41, a42, a43,
        C21, C31, C32, C41, C42, C43, C51, C52, C53, C54,
        gamma, c2, c3, d1, d2, d3,
        h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25)
end

@ROS2(:tableau)

@ROS23(:tableau)

@ROS34PW(:tableau)

@Rosenbrock4(:tableau)

struct RodasTableau{T, T2}
    a::Matrix{T}
    C::Matrix{T}
    gamma::T
    c::Vector{T2}
    d::Vector{T}
    h21::T
    h22::T
    h23::T
    h24::T
    h25::T
    h31::T
    h32::T
    h33::T
    h34::T
    h35::T
end

function Rodas4Tableau(T, T2)
    gamma = convert(T, 1 // 4)
    #BET2P=0.0317D0
    #BET3P=0.0635D0
    #BET4P=0.3438D0
    a = [
        0 convert(T, 1.544000000000000)  0    0     0;
        0 convert(T, 0.9466785280815826) convert(T, 0.2557011698983284)   0     0;
        0 convert(T, 3.314825187068521) convert(T, 2.896124015972201)  convert(T, 0.9986419139977817)    0;
        0 convert(T, 1.221224509226641) convert(T, 6.019134481288629)  convert(T, 12.53708332932087)   -convert(T, 0.6878860361058950)
    ]

    C = [
        0    0     0    0    0;  
        -convert(T, 5.668800000000000)  0     0    0    0;  
        -convert(T, 2.430093356833875)  -convert(T, 0.2063599157091915)   0    0    0;  
        -convert(T, 0.1073529058151375) -convert(T, 9.594562251023355)   -convert(T, 20.47028614809616)  0    0;  
        convert(T, 7.496443313967647) -convert(T, 10.24680431464352)   -convert(T, 33.99990352819905)  convert(T, 11.70890893206160)  0;  
        convert(T, 8.083246795921522)  -convert(T, 7.981132988064893)   -convert(T, 31.52159432874371)  convert(T, 16.31930543123136)  -convert(T, 6.058818238834054)
    ]

    c = [convert(T2, 0.386), convert(T2, 0.21), convert(T2, 0.63)]

    d = [convert(T, 0.2500000000000000), -convert(T, 0.1043000000000000), convert(T, 0.1035000000000000), -convert(T, 0.03620000000000023)]

    h21 = convert(T, 10.12623508344586)
    h22 = -convert(T, 7.487995877610167)
    h23 = -convert(T, 34.80091861555747)
    h24 = -convert(T, 7.992771707568823)
    h25 = convert(T, 1.025137723295662)
    h31 = -convert(T, 0.6762803392801253)
    h32 = convert(T, 6.087714651680015)
    h33 = convert(T, 16.43084320892478)
    h34 = convert(T, 24.76722511418386)
    h35 = -convert(T, 6.594389125716872)

    RodasTableau(a, C, gamma, c, d,
    h21, h22, h23, h24, h25, h31, h32, h33, h34, h35)
end

function Rodas42Tableau(T, T2)
    gamma = convert(T, 1 // 4)
    #BET2P=0.0317D0
    #BET3P=0.0047369D0
    #BET4P=0.3438D0
    a = [
        0     convert(T, 1.402888400000000)  0     0     0;
        0     convert(T, 0.6581212688557198)  -convert(T, 1.320936088384301)  0     0;
        0     convert(T, 7.131197445744498)  convert(T, 16.02964143958207)  -convert(T, 5.561572550509766)   0;
        0     convert(T, 22.73885722420363)  convert(T, 67.38147284535289)  -convert(T, 31.21877493038560)  convert(T, 0.7285641833203814)
    ]

    C = [
        0     -convert(T, 5.104353600000000)  0     0     0         0;  
        0     -convert(T, 2.899967805418783)  convert(T, 4.040399359702244)  0     0          0;  
        0     -convert(T, 32.64449927841361)  -convert(T, 99.35311008728094)  convert(T, 49.99119122405989)   0          0;  
        0     -convert(T, 76.46023087151691)  -convert(T, 278.5942120829058)  convert(T, 153.9294840910643)  convert(T, 10.97101866258358)        0;  
        0     -convert(T, 76.29701586804983)  -convert(T, 294.2795630511232)  convert(T, 162.0029695867566)  convert(T, 23.65166903095270)  -convert(T, 7.652977706771382);
    ]

    c = [convert(T2, 0.3507221), convert(T2, 0.2557041), convert(T2, 0.6817790)]

    d = [convert(T, 0.2500000000000000), -convert(T, 0.06902209999999998), -convert(T, 0.0009671999999999459), -convert(T, 0.08797900000000025)]

    h21 = -convert(T, 38.71940424117216)
    h22 = -convert(T, 135.8025833007622)
    h23 = convert(T, 64.51068857505875)
    h24 = -convert(T, 4.192663174613162)
    h25 = -convert(T, 2.531932050335060)
    h31 = -convert(T, 14.99268484949843)
    h32 = -convert(T, 76.30242396627033)
    h33 = convert(T, 58.65928432851416)
    h34 = convert(T, 16.61359034616402)
    h35 = -convert(T, 0.6758691794084156)

    RodasTableau(a, C, gamma, c, d,
    h21, h22, h23, h24, h25, h31, h32, h33, h34, h35)
end

function Rodas4PTableau(T, T2)
    gamma = convert(T, 1 // 4)
    #BET2P=0.D0
    #BET3P=c3*c3*(c3/6.d0-GAMMA/2.d0)/(GAMMA*GAMMA)
    #BET4P=0.3438D0
    a = [
        0     convert(T, 3)  0     0     0;
        0     convert(T, 1.831036793486759)  convert(T, 0.4955183967433795)  0     0;
        0     convert(T, 2.304376582692669)  -convert(T, 0.05249275245743001)  -convert(T, 1.176798761832782)   0;
        0     -convert(T, 7.170454962423024)  -convert(T, 4.741636671481785)  -convert(T, 16.31002631330971)  -convert(T, 1.062004044111401);
    ]

    C = [
        0     0     0     0     0;  
        -convert(T, 12)   0     0     0     0;  
        -convert(T, 8.791795173947035)   -convert(T, 2.207865586973518)   0     0     0;  
        convert(T, 10.81793056857153)    convert(T, 6.780270611428266)  convert(T, 19.53485944642410)  0     0;  
        convert(T, 34.19095006749676)    convert(T, 15.49671153725963)  convert(T, 54.74760875964130)  convert(T, 14.16005392148534)   0;  
        convert(T, 34.62605830930532)    convert(T, 15.30084976114473)  convert(T, 56.99955578662667)  convert(T, 18.40807009793095)   -convert(T, 5.714285714285717);
    ]

    c = [convert(T2, 3 * gamma), convert(T2, 0.21), convert(T2, 0.63)]

    d = [convert(T, 0.2500000000000000), convert(T, -0.5000000000000000), convert(T, -0.0235040000000000), convert(T, -0.0362000000000000)]

    h21 = convert(T, 25.09876703708589)
    h22 = convert(T, 11.62013104361867)
    h23 = convert(T, 28.49148307714626)
    h24 = -convert(T, 5.664021568594133)
    h25 = convert(T, 0)
    h31 = convert(T, 1.638054557396973)
    h32 = -convert(T, 0.7373619806678748)
    h33 = convert(T, 8.477918219238990)
    h34 = convert(T, 15.99253148779520)
    h35 = -convert(T, 1.882352941176471)

    RodasTableau(a, C, gamma, c, d, 
    h21, h22, h23, h24, h25, h31, h32, h33, h34, h35)
end

function Rodas4P2Tableau(T, T2)
    gamma = convert(T, 1 // 4)

    a = [
        0     convert(T, 3.000000000000000)  0     0     0;
        0     convert(T, 0.906377755268814)  -convert(T, 0.189707390391685)  0     0;
        0     convert(T, 3.758617027739064)  convert(T, 1.161741776019525)  -convert(T, 0.849258085312803)   0;
        0     convert(T, 7.089566927282776)  convert(T, 4.573591406461604)  -convert(T, 8.423496976860259)  -convert(T, 0.959280113459775);
    ] 

    C = [
        0     convert(T, -12.00000000000000)  0     0     0       0;  
        0     convert(T, -6.354581592719008)  convert(T, 0.338972550544623)  0     0        0;  
        0     convert(T, -8.575016317114033)  convert(T, -7.606483992117508)  convert(T, 12.224997650124820)   0          0;  
        0     convert(T, -5.888975457523102)  convert(T, -8.157396617841821)  convert(T, 24.805546872612922)  convert(T, 12.790401512796979)             0;  
        0     convert(T, -4.408651676063871)  convert(T, -6.692003137674639)  convert(T, 24.625568527593117)  convert(T, 16.627521966636085)  convert(T, -5.714285714285718);
    ]

    c = [convert(T2, 0.750000000000000), convert(T2, 0.321448134013046), convert(T2, 0.519745732277726)]

    d = [convert(T, 0.250000000000000), convert(T, -0.500000000000000), convert(T, -0.189532918363016), convert(T, 0.085612108792769)] 

    h21 = convert(T, -5.323528268423303)
    h22 = convert(T, -10.042123754867493)
    h23 = convert(T, 17.175254928256965)
    h24 = convert(T, -5.079931171878093)
    h25 = convert(T, -0.016185991706112)
    h31 = convert(T, 6.984505741529879)
    h32 = convert(T, 6.914061169603662)
    h33 = convert(T, -0.849178943070653)
    h34 = convert(T, 18.104410789349338)
    h35 = convert(T, -3.516963011559032)

    RodasTableau(a, C, gamma, c, d,
    h21, h22, h23, h24, h25, h31, h32, h33, h34, h35)
end

struct Rodas5Tableau{T, T2}
    a::Matrix{T}
    C::Matrix{T}
    gamma::T2
    d::Vector{T}
    c::Vector{T}
    h21::T
    h22::T
    h23::T
    h24::T
    h25::T
    h26::T
    h27::T
    h28::T
    h31::T
    h32::T
    h33::T
    h34::T
    h35::T
    h36::T
    h37::T
    h38::T
    h41::T
    h42::T
    h43::T
    h44::T
    h45::T
    h46::T
    h47::T
    h48::T
end

function Rodas5Tableau(T, T2)
    gamma = convert(T2, 0.19)
    a = [
        0     convert(T, 2.0)  0     0     0     0;
        0     convert(T, 3.040894194418781)  convert(T, 1.041747909077569)  0     0     0;
        0     convert(T, 2.576417536461461)  convert(T, 1.622083060776640)  convert(T, -0.9089668560264532)   0     0;
        0     convert(T, 2.760842080225597)  convert(T, 1.446624659844071)  convert(T, -0.3036980084553738) convert(T, 0.2877498600325443)   0;
        0     convert(T, -14.09640773051259)  convert(T, 6.925207756232704) convert(T, -41.47510893210728)  convert(T, 2.343771018586405)  convert(T, 24.13215229196062);
    ]

    C = [
        0     convert(T, -10.31323885133993)  0     0     0     0     0     0;
        0     convert(T, -21.04823117650003)  convert(T, -7.234992135176716)  0     0     0     0       0;
        0     convert(T, 32.22751541853323)  convert(T, -4.943732386540191)  convert(T, 19.44922031041879)  0     0     0       0;
        0     convert(T, -20.69865579590063)  convert(T, -8.816374604402768)  convert(T, 1.260436877740897)  convert(T, -0.7495647613787146)  0     0       0;
        0     convert(T, -46.22004352711257)  convert(T, -17.49534862857472)  convert(T, -289.6389582892057)  convert(T, 93.60855400400906)  convert(T, 318.3822534212147)  0       0;
        0     convert(T, 34.20013733472935)   convert(T, -14.15535402717690)  convert(T, 57.82335640988400)  convert(T, 25.83362985412365)  convert(T, 1.408950972071624)  convert(T, -6.551835421242162)       0;
        0     convert(T, 42.57076742291101)  convert(T, -13.80770672017997)  convert(T, 93.98938432427124)  convert(T, 18.77919633714503)  convert(T, -31.58359187223370)  convert(T, -6.685968952921985) convert(T, -5.810979938412932);
    ]
    c = [convert(T2, 0.38), convert(T2, 0.3878509998321533), convert(T2, 0.4839718937873840), convert(T2, 0.4570477008819580)]
    d = [convert(T, gamma), convert(T, -0.1823079225333714636), convert(T, -0.319231832186874912), convert(T, 0.3449828624725343), convert(T, -0.377417564392089818)]

    h21 = convert(T, 27.354592673333357)
    h22 = convert(T, -6.925207756232857)
    h23 = convert(T, 26.40037733258859)
    h24 = convert(T, 0.5635230501052979)
    h25 = convert(T, -4.699151156849391)
    h26 = convert(T, -1.6008677469422725)
    h27 = convert(T, -1.5306074446748028)
    h28 = convert(T, -1.3929872940716344)

    h31 = convert(T, 44.19024239501722)
    h32 = convert(T, 1.3677947663381929e-13)
    h33 = convert(T, 202.93261852171622)
    h34 = convert(T, -35.5669339789154)
    h35 = convert(T, -181.91095152160645)
    h36 = convert(T, 3.4116351403665033)
    h37 = convert(T, 2.5793540257308067)
    h38 = convert(T, 2.2435122582734066)

    h41 = convert(T, -44.0988150021747)
    h42 = convert(T, -5.755396159656812e-13)
    h43 = convert(T, -181.26175034586677)
    h44 = convert(T, 56.99302194811676)
    h45 = convert(T, 183.21182741427398)
    h46 = convert(T, -7.480257918273637)
    h47 = convert(T, -5.792426076169686)
    h48 = convert(T, -5.32503859794143)

    # println("---Rodas5---")

    #=
    a71 = -14.09640773051259
    a72 = 6.925207756232704
    a73 = -41.47510893210728
    a74 = 2.343771018586405
    a75 = 24.13215229196062
    a76 = convert(T,1)
    a81 = -14.09640773051259
    a82 = 6.925207756232704
    a83 = -41.47510893210728
    a84 = 2.343771018586405
    a85 = 24.13215229196062
    a86 = convert(T,1)
    a87 = convert(T,1)
    b1 = -14.09640773051259
    b2 = 6.925207756232704
    b3 = -41.47510893210728
    b4 = 2.343771018586405
    b5 = 24.13215229196062
    b6 = convert(T,1)
    b7 = convert(T,1)
    b8 = convert(T,1)
    =#

    RodasTableau(a, C, gamma, d, c,
    h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37,
    h38, h41, h42, h43, h44, h45, h46, h47, h48)
end

function Rodas5PTableau(T, T2)
    gamma = convert(T2, 0.21193756319429014)

    a = [
        0     convert(T, 3.0)  0     0     0            0;
        0     convert(T, 2.849394379747939)  convert(T, 0.45842242204463923)  0     0        0;
        0     convert(T, -6.954028509809101)  convert(T, 2.489845061869568)  convert(T, -10.358996098473584)  0        0;
        0     convert(T, 2.8029986275628964)  convert(T, 0.5072464736228206)  convert(T, -0.3988312541770524)  convert(T, -0.04721187230404641)         0;
        0     convert(T, -7.502846399306121)  convert(T, 2.561846144803919)  convert(T, -11.627539656261098)  convert(T, -0.18268767659942256)  convert(T, 0.030198172008377946);
    ]

    C = [
        0  convert(T, -14.155112264123755)  0   0   0   0   0   0;
        0  convert(T, -17.97296035885952)  convert(T, -2.859693295451294) 0   0   0   0     0;
        0  convert(T, 147.12150275711716)  convert(T, -1.41221402718213) convert(T, 71.68940251302358) 0   0   0    0;
        0  convert(T, 165.43517024871676)  convert(T, -0.4592823456491126) convert(T, 42.90938336958603) convert(T, -5.961986721573306) 0   0    0;
        0  convert(T, 24.854864614690072)  convert(T, -3.0009227002832186) convert(T, 47.4931110020768) convert(T, 5.5814197821558125) convert(T, -0.6610691825249471) 0   0;
        0  convert(T, 30.91273214028599)  convert(T, -3.1208243349937974) convert(T, 77.79954646070892) convert(T, 34.28646028294783) convert(T, -19.097331116725623) convert(T, -28.087943162872662)     0;
        0  convert(T, 37.80277123390563)  convert(T, -3.2571969029072276) convert(T, 112.26918849496327) convert(T, 66.9347231244047) convert(T, -40.06618937091002) convert(T, -54.66780262877968) convert(T, -9.48861652309627);
    ]

    c = [convert(T2, 0.6358126895828704), convert(T2, 0.4095798393397535), convert(T2, 0.9769306725060716), convert(T2, 0.4288403609558664)]

    d = [convert(T, 0.21193756319429014), convert(T, -0.42387512638858027), convert(T, -0.3384627126235924), convert(T, 1.8046452872882734), convert(T, 2.325825639765069)]

    h21 = convert(T, 25.948786856663858)
    h22 = convert(T, -2.5579724845846235)
    h23 = convert(T, 10.433815404888879)
    h24 = convert(T, -2.3679251022685204)
    h25 = convert(T, 0.524948541321073)
    h26 = convert(T, 1.1241088310450404)
    h27 = convert(T, 0.4272876194431874)
    h28 = convert(T, -0.17202221070155493)

    h31 = convert(T, -9.91568850695171)
    h32 = convert(T, -0.9689944594115154)
    h33 = convert(T, 3.0438037242978453)
    h34 = convert(T, -24.495224566215796)
    h35 = convert(T, 20.176138334709044)
    h36 = convert(T, 15.98066361424651)
    h37 = convert(T, -6.789040303419874)
    h38 = convert(T, -6.710236069923372)

    h41 = convert(T, 11.419903575922262)
    h42 = convert(T, 2.8879645146136994)
    h43 = convert(T, 72.92137995996029)
    h44 = convert(T, 80.12511834622643)
    h45 = convert(T, -52.072871366152654)
    h46 = convert(T, -59.78993625266729)
    h47 = convert(T, -0.15582684282751913)
    h48 = convert(T, 4.883087185713722)

    RodasTableau(a,
        C,
        gamma, d, c,
        h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37,
        h38, h41, h42, h43, h44, h45, h46, h47, h48)
end

@RosenbrockW6S4OS(:tableau)

#=
# alpha_ij
A = [0 0 0 0 0 0 0 0
     big"0.38" 0 0 0 0 0 0 0
     big"0.1899188971074152"    big"0.1979321027247381"  0 0 0 0 0 0
     big"0.1110729281178426"    big"0.5456026683145674"  big"-0.1727037026450261" 0 0 0 0 0
     big"0.2329444418850307"    big"0.025099380960713898" big"0.1443314046300300"  big"0.054672473406183418" 0 0 0 0
     big"-0.036201017843430883" big"4.208448872731939"   big"-7.549674427720996"  big"-0.2076823626400282" big"4.585108935472517" 0 0 0
     big"7.585261698003052"     big"-15.57426208319938"  big"-8.814406895608121"  big"1.534698996826085"   big"16.07870828397837" big"0.19" 0 0
     big"0.4646018839086969"    big"0"                   big"-1.720907508837576"  big"0.2910480220957973"  big"1.821778861539924" big"-0.046521258706842056" big"0.19" 0]
# bi
B = [big"0.464601884",0,big"-1.72090751",big"0.29104802",big"1.82177886",big"-0.02674488",big"-0.01977638",big"0.19"]

# Beta_ij

Beta = [big"0.19" 0 0 0 0 0 0 0
        big"0.0076920774666285364"  big"0.19" 0 0 0 0 0 0
        big"-0.058129718999580252"  big"-0.063251113355141360"  big"0.19" 0 0 0 0 0
        big"0.7075715596134048"     big"-0.5980299539145789"   big"0.5294131505610923"    big"0.19" 0 0 0 0
        big"-0.034975026573934865"  big"-0.1928476085817357"    big"0.089839586125126941"  big"0.027613185520411822" big"0.19" 0 0 0
        big"7.585261698003052"      big"-15.57426208319938"     big"-8.814406895608121"    big"1.534698996826085"    big"16.07870828397837"  big"0.19" 0 0
        big"0.4646018839086969"     0                           big"-1.720907508837576"    big"0.2910480220957973"   big"1.821778861539924"  big"-0.046521258706842056" big"0.19" 0
        big"0.4646018839086969"     0                           big"-1.720907508837576"    big"0.2910480220957973"   big"1.821778861539924"  big"-0.026744882930135193" big"-0.019776375776706864" big"0.19"]

Gamma = Beta - A
a = A*Gamma
m = B'*inv(Gamma) # = b_i
C = inv(Diagonal(diag(Gamma))) - inv(Gamma)
c = sum(A,2)
D = Beta*inv(Gamma)
d = sum(Gamma,2)

# Dense output
D_i = tanspose(D)*k_i
y1(Θ) = y₀*(1-Θ) + Θ*(y₁ + (Θ-1)*(D₂ + D₄ + (Θ + 1)*(D₃ + ΘD₄)))

# Determining coefficients
gamma = 0.19
c3 = 0.3878509998321533 == alpha3
c4 = 0.4839718937873840 == alpha4
c5 = 0.4570477008819580 == alpha5
beta3 = 6.8619167645278386e-2
beta4 = 0.8289547562599182
beta5 = 7.9630136489868164e-2
alpha64 = -0.2076823627400282
=#
