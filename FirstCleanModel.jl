# Optimization of a PV - Battery - Heat Pump - Thermal Storage system of an energy efficient residential building.
# This code is an intellectual property of Yuanbei "Fred" Fan, and it is owned by Yuanbei "Fred" Fan, Timothy Hall,
# and Milo Petropolous.

# Yuanbei "Fred" Fan
# CEE Atmosphere and Energy
# Stanford University

######### Initialize tools #########

import Pkg;
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Plots")
Pkg.add("PlotlyJS")
# Initialize JuMP to allow mathematical programming models
using JuMP
using CSV
using DataFrames
using Clp
using PlotlyJS

# Main Code

# Retrieve Data from table
df = CSV.read("C:\\Users\\Fred\\Desktop\\InputTable.csv", DataFrame)

# Set Constant Variable Values
UA = 583; # Asheville House
SetPointT_Low = 68; # [F]
SetPointT_High = 75; # [F]
SetPointW_Low = 0.4; # LEED Regulation
SetPointW_High = 0.6; # LEED Regulation

PeakLightingDensity = 1.705; # [BTU/HR/SF] LA Building 
Area = 2700; # [SF] Ashvielle House
SHGC_diffuse = 0.51; # LA Building 
WindowArea = 338.85; # [SF] Ashvielle House

PeakPlugLoad = 1.705; # [BTU/HR/SF] LA Building 
MaxOccupancy = 4; # [PPL] Ashvielle House

Ventilation = 15; # [CFM/PPL] LA Building
CFMInf = 141.75; # n * V [CFM/PPL] Ashvielle 

PersonLatentHeat = 200; # [BTU/HR/PPL]
PersonSensibleHeat = 300; # [BTU/HR/PPL]
TotalPersonHeat = PersonLatentHeat + PersonSensibleHeat; # [BTU/HR/PPL]

MaxCoolingLoad = 3000; # [BTU/HR/PPL]
MaxHeatingLoad = 3000; # [BTU/HR/PPL]

# Pump Sizing
PumpEfficiency = 0.85; # LA Building
PumpPressureDrop = 50;
PumpReturnTempCooling = 72; # [F] LA Building
PumpSupplyTempCooling = 59; # [F] LA Building
PumpReturnTempHeating = 68; # [F] LA Building
PumpSupplyTempHeating = 104; # [F] LA Building
PumpFlowRateCooling = MaxCoolingLoad / (500 * (PumpReturnTempCooling - PumpSupplyTempCooling));
PumpFlowRateHeating = MaxHeatingLoad / (500 * (PumpReturnTempHeating - PumpSupplyTempHeating));
PumpPowerCooling = PumpFlowRateCooling * PumpPressureDrop / (3960 * PumpEfficiency); # [KW]
PumpPowerHeating = PumpFlowRateHeating * PumpPressureDrop / (3960 * PumpEfficiency); # [KW]

# Fan Sizing
FanEfficiency = 0.65; # LA Building
FanPressureDrop = 3; 
FanReturnTempCooling = 72; # [F] LA Building
FanSupplyTempCooling = 14; # [F] Ashvielle min(TempAmbient)
FanReturnTempHeating = 68; # [F] LA Building
FanSupplyTempHeating = 90; # [F] Ashvielle max(TempAmbient)
FanFlowRateCooling = MaxCoolingLoad / (1.08 * (FanReturnTempCooling - FanSupplyTempCooling));
FanFlowRateHeating = MaxHeatingLoad / (1.08 * (FanReturnTempHeating - FanSupplyTempHeating));
FanPowerCooling = FanFlowRateCooling * FanPressureDrop / (6345 * FanEfficiency); # [KW]
FanPowerHeating = FanFlowRateHeating * FanPressureDrop / (6345 * FanEfficiency); # [KW]
FanSize = max(FanPowerCooling, FanPowerHeating); # [KW]   

ChillerVolume = 80 # [GALLON]
ChillerTemp = 39.2 # [F]
HeaterVolume = 80 # [GALLON]
HeaterTemp = 149 # [F]

ChillSize = (68 - ChillerTemp) * ChillerVolume * 1.001 * 834 # [BTU]
HeatSize = (HeaterTemp - 68) * HeaterVolume * 1.001 * 834 # [BTU]  

BatteryOutLimit = 10 # [KW] change
BatteryInLimit = 10 # [KW] change
BatterySize = 40 # [KWH] change

TotalIterations = 716;

global CurrentStates = zeros(3);

Costs = zeros(TotalIterations);
NoStorageCosts = zeros(TotalIterations);
DeltaCosts = zeros(TotalIterations);
AccumulatedCosts = zeros(TotalIterations);
AccumulatedNoStorageCosts = zeros(TotalIterations);

NPVT = zeros(TotalIterations);
NeHAVCT = zeros(TotalIterations);
NHLT = zeros(TotalIterations);

global TotalCost = 0;
global TotalNoStorageCost = 0;

ActionMatrix = zeros(TotalIterations, 15);


function WindowHeatGain(DAY, LST, DHI, DNI)
    # Given parameter of the location
    LAT = 35.42 # [DEGREES]
    LON = -82.56 # [DEGREES]
    StandardMeridian = -120 # [DEGREES]

    # SHGC vs Incident Angle fitted coefficents [1]
    kA = -0.0000000252
    kB = 0.000002029
    kC = -0.00006709
    kD = 0.0002814
    kE = 0.6

    # Window parameters array (contains different types in each direction)
    WindowQuantityN = [4, 2]
    WindowHeightN = [4, 3]
    WindowWidthN = [4, 4]
    OverhangN = [1, 1] # change

    WindowQuantityE = [0]
    WindowHeightE = 0
    WindowWidthE = 0
    OverhangE = 0

    WindowQuantityS = [8, 1, 1, 1] # 1 sunspace accounted at the end
    WindowHeightS = [6, 6, 4, 8]
    WindowWidthS = [40/12, 9, 4, 18]
    OverhangS = [1, 1, 1, 1] # change

    WindowQuantityW = [1, 1]
    WindowHeightW = [4, 4]
    WindowWidthW = [3.046875, 26/12]
    OverhangW = [1, 1] # change

    # Calculate the Window Area from each direction
    global TotalWindowAreaN = 0
    for i in eachindex(WindowQuantityN)
        WindowAreaN = WindowQuantityN[i] * WindowHeightN[i] * WindowWidthN[i]   
        global TotalWindowAreaN += WindowAreaN
    end

    global TotalWindowAreaE = 0
    for i in eachindex(WindowQuantityE)
        WindowAreaE = WindowQuantityE[i] * WindowHeightE[i] * WindowWidthE[i]   
        global TotalWindowAreaE += WindowAreaE
    end

    global TotalWindowAreaS = 0
    for i in eachindex(WindowQuantityS)
        WindowAreaS = WindowQuantityS[i] * WindowHeightS[i] * WindowWidthS[i]   
        global TotalWindowAreaS += WindowAreaS
    end

    global TotalWindowAreaW = 0
    for i in eachindex(WindowQuantityW)
        WindowAreaW = WindowQuantityW[i] * WindowHeightW[i] * WindowWidthW[i]   
        global TotalWindowAreaW += WindowAreaW
    end

    Γ = deg2rad((360 / 365) * (DAY - 1)) # [RAD]

    # Calculate the Equation of Time (ET) [MIN]
    EQT = 2.2918 * (0.0075 + 0.1868 * cos(Γ) - 3.2077 * sin(Γ) - 1.4615 * cos(2 * Γ) - 4.089 * sin(2 * Γ))

    # Calculate the Apparent Solar Time (AST) [HR]
    ApparentTime = LST + (EQT / 60) + ((LON - StandardMeridian) / 15) # LST is the local standard hour

    # Calculate the Hour Angle (H) [DEGREES]
    HourAngle = 15 * (ApparentTime - 12)

    # Calculate the Declination (σ) [DEGREES]
    Declination = 23.45 * sin(deg2rad((360 / 365) * (DAY + 284)))

    # Calculate the Solar Altitude (β) [DEGREES]
    SolarAltitude = rad2deg(asin(cos(deg2rad(LAT)) * cos(deg2rad(Declination)) * cos(deg2rad(HourAngle)) + sin(deg2rad(LAT)) * sin(deg2rad(Declination))))

    # Calculate the Beam Normal Irradiance (Eb) [BTU/HR/SF]
    Eb = DNI / pi
    # Calculate the Diffuse Horizontal Irradiance (Ed) [BTU/HR/SF]
    Ed = DHI / pi

    # Calculate the Solar Azimuth (Φ) [DEGREES]
    if ApparentTime < 12
        SolarAz = rad2deg(acos((sin(deg2rad(SolarAltitude))*sin(deg2rad(LAT))-sin(deg2rad(Declination)))/(cos(deg2rad(SolarAltitude))*cos(deg2rad(LAT)))))*-1
    else
        SolarAz = rad2deg(acos((sin(deg2rad(SolarAltitude))*sin(deg2rad(LAT))-sin(deg2rad(Declination)))/(cos(deg2rad(SolarAltitude))*cos(deg2rad(LAT)))))
    end

    # Calculate the Surface Solar Azimuth (γ) for all 4 directions [DEGREES]
    SurfaceSolarAzN = SolarAz - 180
    SurfaceSolarAzE = SolarAz + 90
    SurfaceSolarAzS = SolarAz
    SurfaceSolarAzW = SolarAz - 90

    # Calculate the Incident Angle (θ) for all 4 directions [DEGREES]
    IAN = rad2deg(acos(cos(deg2rad(SolarAltitude)) * cos(deg2rad(SurfaceSolarAzN))))
    IAE = rad2deg(acos(cos(deg2rad(SolarAltitude)) * cos(deg2rad(SurfaceSolarAzE))))
    IAS = rad2deg(acos(cos(deg2rad(SolarAltitude)) * cos(deg2rad(SurfaceSolarAzS))))
    IAW = rad2deg(acos(cos(deg2rad(SolarAltitude)) * cos(deg2rad(SurfaceSolarAzW))))

    # Calculate the Vertical to Horizontal Irradiance Ratio (Y) for all 4 directions [1]
    V2HIRN = max(0.45, 0.55 + 0.437 * cos(deg2rad(IAN)) + 0.313 * cos(deg2rad(IAN))^2)
    V2HIRE = max(0.45, 0.55 + 0.437 * cos(deg2rad(IAE)) + 0.313 * cos(deg2rad(IAE))^2)
    V2HIRS = max(0.45, 0.55 + 0.437 * cos(deg2rad(IAS)) + 0.313 * cos(deg2rad(IAS))^2)
    V2HIRW = max(0.45, 0.55 + 0.437 * cos(deg2rad(IAW)) + 0.313 * cos(deg2rad(IAW))^2)

    # Calculate the Profile Angle (Ω) for all 4 directions [DEGREES]
    ProfileAngleN = rad2deg(atan(tan(deg2rad(SolarAltitude))/cos(deg2rad(SurfaceSolarAzN))))
    ProfileAngleE = rad2deg(atan(tan(deg2rad(SolarAltitude))/cos(deg2rad(SurfaceSolarAzE))))
    ProfileAngleS = rad2deg(atan(tan(deg2rad(SolarAltitude))/cos(deg2rad(SurfaceSolarAzS))))
    ProfileAngleW = rad2deg(atan(tan(deg2rad(SolarAltitude))/cos(deg2rad(SurfaceSolarAzW))))

    # Calculate the Beam Component Irradiance (Etb) for all 4 directions [BTU/HR/SF]
    EtbN = 0;
    if SurfaceSolarAzN >= -90 && SurfaceSolarAzN <= 90 && SolarAltitude >= 0
        EtbN = cos(deg2rad(ProfileAngleN)) * Eb
    else
        EtbN = 0
    end
    
    EtbE = 0;
    if SurfaceSolarAzE >= -90 && SurfaceSolarAzE <= 90 && SolarAltitude >= 0
        EtbE = cos(deg2rad(ProfileAngleE)) * Eb
    else
        EtbE = 0
    end
    
    EtbS = 0;
    if SurfaceSolarAzS >= -90 && SurfaceSolarAzS <= 90 && SolarAltitude >= 0
        EtbS = cos(deg2rad(ProfileAngleS)) * Eb
    else
        EtbS = 0
    end
    
    EtbW = 0;
    if SurfaceSolarAzW >= -90 && SurfaceSolarAzW <= 90 && SolarAltitude >= 0
        EtbW = cos(deg2rad(ProfileAngleW)) * Eb
    else
        EtbW = 0
    end

    # Calculate the SHGC(Adjusted for Sun Angle) SHGC(θ) for all 4 directions [1]
    SHGCN =  max(kA * IAN ^ 4 + kB * IAN ^ 3 + kC * IAN ^ 2 + kD * IAN  + kE, 0)
    SHGCE =  max(kA * IAE ^ 4 + kB * IAE ^ 3 + kC * IAE ^ 2 + kD * IAE  + kE, 0)
    SHGCS =  max(kA * IAS ^ 4 + kB * IAS ^ 3 + kC * IAS ^ 2 + kD * IAS  + kE, 0)
    SHGCW =  max(kA * IAW ^ 4 + kB * IAW ^ 3 + kC * IAW ^ 2 + kD * IAW  + kE, 0)

    # Calculate the Window Sunlit Areas for all 4 directions [SF]
    global TotalSunlitAreaN = 0
    for i in eachindex(WindowQuantityN)
        SunlitAreaN = 0
        if SurfaceSolarAzN >= -90 && SurfaceSolarAzN <= 90 && SolarAltitude >= 0 && ((WindowHeightN[i] - (OverhangN[i] * tan(deg2rad(ProfileAngleN)))) * WindowWidthN[i]) > 0
            SunlitAreaN = (WindowHeightN[i] - (OverhangN[i] * tan(deg2rad(ProfileAngleN)))) * WindowWidthN[i] * WindowQuantityN[i]
        end
        global TotalSunlitAreaN += SunlitAreaN
    end

    global TotalSunlitAreaE = 0
    for i in eachindex(WindowQuantityE)
        SunlitAreaE = 0
        if SurfaceSolarAzE >= -90 && SurfaceSolarAzE <= 90 && SolarAltitude >= 0 && ((WindowHeightE[i] - (OverhangE[i] * tan(deg2rad(ProfileAngleE)))) * WindowWidthE[i]) > 0
            SunlitAreaE = (WindowHeightE[i] - (OverhangE[i] * tan(deg2rad(ProfileAngleE)))) * WindowWidthE[i] * WindowQuantityE[i]
        end
        global TotalSunlitAreaE += SunlitAreaE
    end

    global TotalSunlitAreaS = 0
    for i in eachindex(WindowQuantityS)
        SunlitAreaS = 0
        if SurfaceSolarAzS >= -90 && SurfaceSolarAzS <= 90 && SolarAltitude >= 0 && ((WindowHeightS[i] - (OverhangS[i] * tan(deg2rad(ProfileAngleS)))) * WindowWidthS[i]) > 0
            SunlitAreaS = (WindowHeightS[i] - (OverhangS[i] * tan(deg2rad(ProfileAngleS)))) * WindowWidthS[i] * WindowQuantityS[i]
        end
        global TotalSunlitAreaS += SunlitAreaS
    end

    global TotalSunlitAreaW = 0
    for i in eachindex(WindowQuantityW)
        SunlitAreaW = 0
        if SurfaceSolarAzW >= -90 && SurfaceSolarAzW <= 90 && SolarAltitude >= 0 && ((WindowHeightW[i] - (OverhangW[i] * tan(deg2rad(ProfileAngleW)))) * WindowWidthW[i]) > 0
            SunlitAreaW = (WindowHeightW[i] - (OverhangW[i] * tan(deg2rad(ProfileAngleW)))) * WindowWidthW[i] * WindowQuantityW[i]
        end
        global TotalSunlitAreaW += SunlitAreaW
    end

    # Calculate the Diffuse Component Irradiance (Etd) for all 4 directions [BTU/HR/SF]
    EtdN = Ed * V2HIRN
    EtdE = Ed * V2HIRE
    EtdS = Ed * V2HIRS
    EtdW = Ed * V2HIRW

    # Calculate the Ground Component Irradiance (Etr) for all 4 directions [BTU/HR/SF]
    EtrN = (Eb * sin(SolarAltitude) + Ed) * 0.2
    EtrE = max((Eb * sin(SolarAltitude) + Ed) * 0.1, 0)
    EtrS = (Eb * sin(SolarAltitude) + Ed) * 0.2
    EtrW = (Eb * sin(SolarAltitude) + Ed) * 0.2

    # Calculate the Windows Direct Solar Heat Gain (qb) for all 4 directions [BTU/HR]
    qbN = EtbN * TotalSunlitAreaN * SHGCN
    qbE = EtbE * TotalSunlitAreaE * SHGCE
    qbS = EtbS * TotalSunlitAreaS * SHGCS
    qbW = EtbW * TotalSunlitAreaW * SHGCW

    # Calculate the Windows Diffuse Solar Heat Gain (qd) for all 4 directions [BTU/HR]
    qdN = (EtdN + EtrN) * SHGC_diffuse * TotalWindowAreaN
    qdE = (EtdE + EtrE) * SHGC_diffuse * TotalWindowAreaE
    qdS = (EtdS + EtrS) * SHGC_diffuse * TotalWindowAreaS
    qdW = (EtdW + EtrW) * SHGC_diffuse * TotalWindowAreaW

    # Calculate the Total Windows Heat Gain (qWindow) [BTU/HR]
    qWindow = qbN + qdN + qbE + qdE + qbS + qdS + qbW + qdW

    return qWindow
end


function GetData(iteration, DataMatrix)
    # Collect data for 168 iterations from the big data table. Return outputs needed for the optimize funciton.
    Results = zeros(168, 6);
    for i = 1:168
        DAY = ceil((iteration + i - 1)/24);
        println(DAY)
        Results[i, 1] = DataMatrix[i, 2]/100 ; # Time (LST) [HR]
        Results[i, 2] = (DataMatrix[i, 3] * 1.8) + 32 ; # Ambient temperature [F]
        
        # Calculate delta T for each row
        TempAmbient = Results[i,2];
        if TempAmbient > SetPointT_High
            TempDelta = TempAmbient - SetPointT_High;
        elseif TempAmbient < SetPointT_Low
            TempDelta = TempAmbient - SetPointT_Low;  # Negative delta T means negative q, which is heating load
        else 
            TempDelta = 0;
        end
        # Calculate qEvenlope
        qEvenlope = UA * TempDelta;
        
        # Heat gain from lighting [BTU/HR]
        PercentLighting = DataMatrix[i, 5];
        qLighting = PeakLightingDensity * Area * PercentLighting;

        # Electricity consumption from lighting [KWH]
        eLighting = qLighting / 3412.14; # * 1[HR]

        # Heat gain from plugs [BTU/HR]
        PercentPlug = DataMatrix[i, 6];
        qPlugs = PeakPlugLoad * Area * PercentPlug;

        # Electricity consumption from plugs [KWH]
        ePlugs = qPlugs / 3412.14; # * 1(HR)

        # Heat gain from occupancy [BTU/HR]
        PercentOccupied = DataMatrix[i, 7];
        qOccupancy = PercentOccupied * MaxOccupancy * TotalPersonHeat;
        
        HumidityAmbient = DataMatrix[i, 4]/100; # Ambient Humidity (RH)
        # Calculate delta W for each row
        if HumidityAmbient > SetPointW_High
            MoistDelta = HumidityAmbient - SetPointW_High;
        elseif HumidityAmbient < SetPointW_Low
            MoistDelta = HumidityAmbient - SetPointW_Low;  
        else 
            MoistDelta = 0;
        end
        # Heat transfer from infiltration + Ventilation (Sensible(temperature difference) + Latent(moisture difference)) [BTU/HR]
        qInfiltration = 1.08 * TempDelta * CFMInf + 4840 * MoistDelta * CFMInf;

        CFMVen = min(5 * PercentOccupied * MaxOccupancy + 0.06 * Area, PercentOccupied * MaxOccupancy * Ventilation); # LA Building
        qVentilation = 1.08 * TempDelta * CFMVen + 4840 * MoistDelta * CFMVen;

        # Heat gain from windows

        DHI = DataMatrix[i, 9];
        DNI = DataMatrix[i, 10];

        qWindow = WindowHeatGain(DAY, DataMatrix[i, 2], DHI, DNI);
        
        # Total load
        qTotal = qEvenlope + qLighting + qPlugs + qOccupancy + qInfiltration + qVentilation + qWindow;
        Results[i, 3] = qTotal;
        
        # Electricity consumption from pump (Distribution of chilled or hot water in radiant system) [KWH]
        ePump = 0;

        # check on concepts
        
        if qTotal > 0 # needs cooling
            ePump = PumpPowerCooling * (qTotal / MaxCoolingLoad) ^ 3; # * 1[HR]
        else  # needs heating 
            ePump = PumpPowerHeating * (qTotal / MaxHeatingLoad) ^ 3; # * 1[HR]
        end
        # Electricity consumption from fan (Distribution of air, DOAS fan energy) [KWH]
        CFMFan = 0;
        if TempDelta == 0
            CFMFan = 0;
        else 
            CFMFan = qVentilation / (1.08 * TempDelta);
        end
        eFan = FanSize * (CFMFan / (max(FanFlowRateCooling, FanFlowRateHeating))) ^ 3; # * 1[HR]
        
        # Total electricity consumption from house [KWH]
        eHouse = eLighting + ePlugs; #+ ePump + eFan;
        Results[i, 4] = eHouse;

        Results[i, 5] = DataMatrix[i, 8]; # TOU Rate
        Results[i, 6] = 2*DataMatrix[i, 11]/1000; # PV Output [KW]
    end
    return Results
end

function Optimize(OptInput, InitialStates)
    #Set timesteps 

    TimeStart = 1;
    TimeEnd = 168;
    TIME = collect(TimeStart:1:TimeEnd); # Collect time steps into a vector
    NumTime = length(TIME); # Number of time steps, useful for indexing

    InitialStorageBattery = InitialStates[1];
    InitialStorageHeat = InitialStates[2];
    InitialStorageChill = InitialStates[3];
       
    TempAmbient = OptInput[:, 2];
    HouseLoad = OptInput[:, 4];
    TOU = OptInput[:, 5];
    PVGeneration = OptInput[:, 6];
      
    HeatingLoad = zeros(NumTime);
    CoolingLoad = zeros(NumTime);  

    for i = 1:168
        if OptInput[i, 3] > 0
            CoolingLoad[i] = OptInput[i, 3];
        else
            HeatingLoad[i] = -OptInput[i, 3];
        # Heating and cooling load are now both defined as positive [BTU/HR]
        end
    end

    ########## Declare model  ##########

    # Define the model name and solver. In this case, model name is "m"
    m = Model(Clp.Optimizer)

    ######## Decision variables ########

    @variable(m, PV2H[1:NumTime] >= 0); # [KWH]

    @variable(m, PV2G[1:NumTime] >= 0); # [KWH]

    @variable(m, PV2B[1:NumTime] >= 0); # [KWH]

    @variable(m, B2G[1:NumTime] >= 0); # [KWH]

    @variable(m, B2H[1:NumTime] >= 0); # [KWH]

    @variable(m, G2B[1:NumTime] >= 0); # [KWH]

    @variable(m, G2H[1:NumTime] >= 0); # [KWH]

    @variable(m, H2HP[1:NumTime] >= 0); # [KWH]

    @variable(m, HP2H[1:NumTime] >= 0); # [BTU]

    @variable(m, TS2H[1:NumTime] >= 0); # [BTU]

    @variable(m, HP2TS[1:NumTime] >= 0); # [BTU]

    @variable(m, H2C[1:NumTime] >= 0); # [KWH]

    @variable(m, C2H[1:NumTime] >= 0); # [BTU]

    @variable(m, CS2H[1:NumTime] >= 0); # [BTU]

    @variable(m, C2CS[1:NumTime] >= 0); # [BTU]]

    @variable(m, InStorageBattery[1:NumTime+1] >= 0); # [KWH]

    @variable(m, InStorageHeat[1:NumTime+1] >= 0); # [BTU]

    @variable(m, InStorageChill[1:NumTime+1] >= 0); # [BTU]
    
    ######## Objective Functions #########
    
    # Single objective for minimizing cost

    @objective(m, Min, sum(TOU[t]*(G2H[t] + G2B[t] - PV2G[t] - B2G[t]) for t=1:NumTime));

    ############# Constraints ############

    # Battery storage initialization constraint
    @constraint(m, InStorageBattery[1] == InitialStorageBattery);

    # Hot water storage initialization constraint
    @constraint(m, InStorageHeat[1] == InitialStorageHeat);

    # Cold water storage initialization constraint
    @constraint(m, InStorageChill[1] == InitialStorageChill);

    # PV power balance constraint, node at PV
    for t = 1:NumTime
        @constraint(m, PVGeneration[t] ==  PV2B[t] + PV2H[t] + PV2G[t]);
    end

    # House electricity load constraint, node at house
    for t = 1:NumTime
        @constraint(m, HouseLoad[t] + H2HP[t] + H2C[t] <= PV2H[t] + B2H[t] + G2H[t]);
    end

    # House heating load constraint, node at house
    for t = 1:NumTime
        @constraint(m, HeatingLoad[t] <= HP2H[t] + TS2H[t]);
    end

    # House cooling load constraint, node at house
    for t = 1:NumTime
        @constraint(m, CoolingLoad[t] <= C2H[t] + CS2H[t]);
    end

    # Battery storage balance constraint, node at battery
    for t = 1:(NumTime - 1)
        @constraint(m, InStorageBattery[t+1] == InStorageBattery[t] + PV2B[t] + G2B[t] - B2G[t] - B2H[t]);
    end

    # Battery Discharging constraint, node at battery
    for t = 1:NumTime
        @constraint(m, B2G[t] + B2H[t] <= InStorageBattery[t]);
    end

    # Battery Power Discharging constraint, node at battery
    for t = 1:NumTime
        @constraint(m, B2G[t] + B2H[t] <= BatteryOutLimit);  #model degradation
    end

    # Battery Power Charging constraint, node at battery
    for t = 1:NumTime
        @constraint(m, PV2B[t] + G2B[t] <= BatteryInLimit);  #model degradation
    end

    # Battery Storage size constraint, node at battery
    for t = 1:NumTime
        @constraint(m, InStorageBattery[t] <= BatterySize); #model degradation
    end

    # Hot water storage balance constraint, node at hot water storage [BTU]
    for t = 1:(NumTime - 1)
        @constraint(m, InStorageHeat[t+1] == InStorageHeat[t] + HP2TS[t] - TS2H[t]); #think about loss
    end

    # Hot water Storage size constraint, node at hot water storage [BTU]
    for t = 1:NumTime
        @constraint(m, InStorageHeat[t] <= HeatSize);
    end

    # Hot water storage availablility constraint, node at hot water storage [BTU]
    for t = 1:NumTime
        @constraint(m, TS2H[t] <= InStorageHeat[t]);
    end
    
    # Cold water storage balance constraint, node at cold water storage [BTU]
    for t = 1:(NumTime - 1)
        @constraint(m, InStorageChill[t+1] == InStorageChill[t] + C2CS[t] - CS2H[t]); #think about loss
    end

    # Cold water Storage size constraint, node at cold water storage [BTU]
    for t = 1:NumTime
        @constraint(m, InStorageChill[t] <= ChillSize);
    end

    # Cold water storage availablility constraint, node at cold water storage [BTU]
    for t = 1:NumTime
        @constraint(m, CS2H[t] <= InStorageChill[t]);
    end
    
    # Air Source Heat Pump energy balance constraint, node at ASHP [KWH to BTU]
    # move numerical values to top of the code
    for t = 1:NumTime
        @constraint(m, H2HP[t] * (3412.14 * (0.0378 * TempAmbient[t] + 1.18557)) ==  HP2H[t] + HP2TS[t]);
    end

    # Chiller energy balance constraint, node at chiller [KWH to BTU]
    for t = 1:NumTime
        @constraint(m, H2C[t] * (3412.14 * (0.0334 * TempAmbient[t] + 5.8451)) ==  C2H[t] + C2CS[t]);
    end

    ########### Print and solve ##########

    # print(m);

    optimize!(m);
    
    # Return all decisions at first time step
    BestActions = zeros(15);
    BestActions[1] = value.(PV2H[1]);
    BestActions[2] = value.(PV2G[1]);
    BestActions[3] = value.(PV2B[1]);
    BestActions[4] = value.(B2G[1]);
    BestActions[5] = value.(B2H[1]);
    BestActions[6] = value.(G2B[1]);
    BestActions[7] = value.(G2H[1]);
    BestActions[8] = value.(H2HP[1]);
    BestActions[9] = value.(HP2H[1]);
    BestActions[10] = value.(TS2H[1]);
    BestActions[11] = value.(HP2TS[1]);
    BestActions[12] = value.(H2C[1]);
    BestActions[13] = value.(C2H[1]);
    BestActions[14] = value.(CS2H[1]);
    BestActions[15] = value.(C2CS[1]);
  
    NetCost = (BestActions[6] + BestActions[7] - BestActions[2] - BestActions[4]) * TOU[1];
    
    CurrentStates = zeros(3);
    CurrentStates[1] = value.(InStorageBattery[2]);
    CurrentStates[2] = value.(InStorageHeat[2]);
    CurrentStates[3] = value.(InStorageChill[2]);
    
    return BestActions, NetCost, CurrentStates
end

function NoStorage(NormalInput)
    # This funciton describes a senario where there is no thermal storage tank, chilled water tank, or battery.
    TempAmbient = NormalInput[2];
    HouseLoad =  NormalInput[4];
    TOU =  NormalInput[5];
    PVGeneration =  NormalInput[6];
    eHVAC = 0;

    if  NormalInput[3] > 0
        CoolingLoad = NormalInput[3];
        eHVAC = CoolingLoad / (3412.14 * (0.0334 * TempAmbient + 5.8451));
    else
        HeatingLoad = -NormalInput[3];
        eHVAC = HeatingLoad / (3412.14 * (0.0378 * TempAmbient + 1.18557));
    end
    # Calculate cost
    Cost = TOU * (HouseLoad + eHVAC - PVGeneration);
    
    return Cost, HouseLoad, eHVAC, PVGeneration
end

for i = 1:TotalIterations
    inputtable = GetData(i, df[i:i+167,:]);
    # Get data for non storage function.
    CurrentRow = inputtable[1, :]
    GoodActions, NetCost, NewStates = Optimize(inputtable, CurrentStates);
    # Get cost for the no storage situation.
    NoStorageCost,nHL, neHVAC, nPV = NoStorage(CurrentRow);
    # Update the current states of heat, chill, and battery storage.
    for x = 1:3
        global CurrentStates[x] = NewStates[x];
    end
    # Record the actions from each iteration
    for j = 1:15
        ActionMatrix[i, j] = GoodActions[j];
    end
    # Add the cost
    Costs[i] = NetCost;
    global TotalCost += NetCost;
    # Compute the cost difference between optimized system and no storage system
    DeltaCost = NoStorageCost - NetCost;
    global TotalNoStorageCost += NoStorageCost;
    NoStorageCosts[i] = NoStorageCost;
    NeHAVCT[i] = neHVAC;
    NPVT[i] = nPV;
    NHLT[i] = nHL;
    DeltaCosts[i] = DeltaCost;
    AccumulatedCosts[i] = sum(Costs[n] for n = 1:i)
    AccumulatedNoStorageCosts[i] = sum(NoStorageCosts[n] for n = 1:i)
end


trace1 = scatter(
    x = 1:TotalIterations,  
    y=AccumulatedCosts,
    name="Optimized Storage System"
)
trace2 = scatter(
    x = 1:TotalIterations, 
    y=AccumulatedNoStorageCosts,
    name="No Storage Senario"
)

#=
trace1 = scatter(
    x = 1:TotalIterations,  
    y=NPVT,
    name="PV"
)

trace2 = scatter(
    x = 1:TotalIterations, 
    y=NHLT,
    name="House Load"
)

trace3 = scatter(
    x = 1:TotalIterations, 
    y=NeHAVCT,
    name="Heat Pump"
)
=#
plot([trace1, trace2], Layout(legend_title_text="Accumulated Energy Cost over time"))




