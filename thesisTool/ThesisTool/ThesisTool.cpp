/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include <ThesisTool/VehicleAerodynamicCoefficients.h>

#include <ThesisTool/applicationOutput.h>

#include <ThesisTool/recoverymodels.h>

Eigen::Vector3d my_func(int t)
{
    return Eigen::Vector3d (1,0,0);
}

Eigen::Vector3d my_func2(void)
{
    return Eigen::Vector3d (1,0,0);
}

int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::root_finders;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::propagators;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;
    using namespace tudat::propulsion;
    using namespace tudat;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3100.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 0.5;

    // Define simulation body settings.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "First Stage" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Upper Stage" ] = std::make_shared< simulation_setup::Body >( );

    // Define inputs
    std::string atmosphereFile = "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";

    std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables;
    dependentVariables.push_back( aerodynamics::pressure_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::temperature_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::specific_heat_ratio_dependent_atmosphere );
    dependentVariables.push_back( aerodynamics::density_dependent_atmosphere );

    double specificGasConstant = 197.0;

    double ratioOfSpecificHeats = 1.3;

    interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_default_value_with_warning;

    double defaultExtrapolationValue = TUDAT_NAN;

    // Set tabulated atmosphere
    bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereFile, dependentVariables, specificGasConstant, ratioOfSpecificHeats, boundaryHandling, defaultExtrapolationValue );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //// USER INPUTS ////
    BodyProperties FirstStage;
    BodyProperties Vehicle;
    BodyProperties UpperStage;

    bool OrbitalLauncher=false;

    InitialProperties InitialConditions;

    InitialConditions.H0=0;
    InitialConditions.V0=1;
    InitialConditions.heading=0;
    InitialConditions.flightpath=85;
    InitialConditions.latitude=0;
    InitialConditions.longitude=0;

    // Create vehicle aerodynamic coefficients
    Vehicle.referenceArea = 0.0607;
    FirstStage.referenceArea = 0.0607;
    UpperStage.referenceArea = 0.0607;

    Vehicle.InitialWetMass=324;
    FirstStage.InitialWetMass=324-16;
    UpperStage.InitialWetMass=16;

    Vehicle.InitialEmptyMass=100;
    FirstStage.InitialEmptyMass=100;
    UpperStage.InitialEmptyMass=16;

    Vehicle.constantAOA=0;
    FirstStage.constantAOA=0;
    UpperStage.constantAOA=0;

    // Set Which recovery System is used This can either be "Propulsive" "Parachute" "Rotor" "Winged"
//    FirstStage.RecoverySystem={"Propulsive","Parachute"};

    Vehicle.CDfile="C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CD.txt";
    Vehicle.CLfile="C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CL.txt";

    FirstStage.CDfile="C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CD.txt";
    FirstStage.CLfile="C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CL.txt";

    UpperStage.CDfile="C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CD.txt";
    UpperStage.CLfile="C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CL.txt";

    std::string launchThrustFile="C:/tudatBundle/thesisTool/ThesisTool/ThrustValues.txt";

    ///////////// SETTINGS FOR AERODYNAMICS VALUES ////////////
    // Define physical meaning of independent variables, in this case Mach number and angle of attack
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );

    // Define reference frame in which the loaded coefficients are defined.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;

    ///////////// DRAG FOR ENTIRE VEHICLE /////////////////////
    // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles1;
    forceCoefficientFiles1[ 0 ] = Vehicle.CDfile; // Set drag coefficient file
    forceCoefficientFiles1[ 2 ] = Vehicle.CLfile; // Set lift coefficient file

    // Load and parse files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings1 =
        readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles1, Vehicle.referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,areCoefficientsInNegativeAxisDirection );

    // Create and set aerodynamic coefficients
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings1, "Vehicle" ) );
    std::cout << "Loaded Aerodynamic Data: Vehicle"<<std::endl;

    ///////////// DRAG FOR FIRST STAGE ///////////////////////
    // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles2;
    forceCoefficientFiles2[ 0 ] = FirstStage.CDfile; // Set drag coefficient file
    forceCoefficientFiles2[ 2 ] = FirstStage.CLfile; // Set lift coefficient file

    // Load and parse files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings2 =
        readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles2, FirstStage.referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,areCoefficientsInNegativeAxisDirection );

    // Create and set aerodynamic coefficients
    bodyMap[ "First Stage" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings2, "First Stage" ) );
    std::cout << "Loaded Aerodynamic Data : First Stage"<<std::endl;
    ///////////// DRAG FOR FIRST STAGE /////////////////////////
    // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles3;
    forceCoefficientFiles3[ 0 ] = UpperStage.CDfile; // Set drag coefficient file
    forceCoefficientFiles3[ 2 ] = UpperStage.CLfile; // Set lift coefficient file

    // Load and parse files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings3 =
        readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles3, UpperStage.referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,areCoefficientsInNegativeAxisDirection );

    // Create and set aerodynamic coefficients
    bodyMap[ "Upper Stage" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings3, "Upper Stage" ) );
    std::cout << "Loaded Aerodynamic Data : Upper Stage"<<std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////        Recovery Models     ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FirstStage.setFixedRecMass();
    UpperStage.FixedRecMass=0;
    /// Creating Total Masses
    FirstStage.setTotalMass();
//    FirstStage.setRecoveryMass();
//    FirstStage.TotalMass=FirstStage.InitialWetMass+FirstStage.RecoveryMass+FirstStage.FixedRecMass;
    if (OrbitalLauncher==false){
        UpperStage.TotalMass=UpperStage.InitialWetMass;
    };
    Vehicle.TotalMass=FirstStage.TotalMass+UpperStage.TotalMass;

    bodyMap[ "Vehicle" ]->setConstantBodyMass( Vehicle.TotalMass );
    bodyMap[ "First Stage" ]->setConstantBodyMass( FirstStage.TotalMass );
    bodyMap[ "Upper Stage" ]->setConstantBodyMass( UpperStage.TotalMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define data to be used for thrust as a function of time.
    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    std::shared_ptr< FromFileDataMapSettings< Eigen::Vector3d > > thrustDataSettings =
            std::make_shared< FromFileDataMapSettings< Eigen::Vector3d > >( launchThrustFile );

    // Direction Based Guidance
    std::function< Eigen::Vector3d(const double) > thrustunitfunction=my_func;
    std::function< Eigen::Vector3d() > bodyunitfunc=my_func2;
    std::make_shared<DirectionBasedForceGuidance> ( thrustunitfunction,"Vehicle",bodyunitfunc);

    // Define interpolator settings.
    std::shared_ptr< InterpolatorSettings > thrustInterpolatorSettings =
            std::make_shared< InterpolatorSettings >( linear_interpolator );

    // Create data interpolation settings
    std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolatorSettings =
            std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                thrustDataSettings, thrustInterpolatorSettings );

    // Define specific impulse
    double constantSpecificImpulse = 210.0;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       thrustDataInterpolatorSettings, constantSpecificImpulse,
                                                       lvlh_thrust_frame, "Earth" ) );
    std::cout << "Thrust Loaded"<<std::endl;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );


    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        bodyMap.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                [=]( ){ return Vehicle.constantAOA; } );
    std::cout << "Acceleration Map Set"<<std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for the vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.
    // Set Initial location of the vehicle
    Eigen::Vector6d InitialLaunchState;
    InitialLaunchState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + InitialConditions.H0;
    InitialLaunchState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            unit_conversions::convertDegreesToRadians( InitialConditions.latitude );
    InitialLaunchState( SphericalOrbitalStateElementIndices::longitudeIndex ) =
            unit_conversions::convertDegreesToRadians( InitialConditions.longitude );
    InitialLaunchState( SphericalOrbitalStateElementIndices::speedIndex ) = InitialConditions.V0;
    InitialLaunchState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            unit_conversions::convertDegreesToRadians( InitialConditions.flightpath );
    InitialLaunchState( SphericalOrbitalStateElementIndices::headingAngleIndex ) =
            unit_conversions::convertDegreesToRadians( InitialConditions.heading );
    // Convert apollo state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                InitialLaunchState );
    std::cout << "Initial Conditions Set"<<std::endl;

    // Convert the state to the global (inertial) frame.
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( "Earth" )->getRotationalEphemeris( );
    systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

    // Define propagation termination conditions (when altitude is 0m).
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable  =
          std::make_shared< SingleDependentVariableSaveSettings >(
              altitude_dependent_variable, "Vehicle", "Earth" );

    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared < propagators::PropagationDependentVariableTerminationSettings>
            (terminationDependentVariable,-1e3,true);


    // Define settings for propagation of translational dynamics.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );

    // Create mass rate models
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = createMassRateModel(
                "Vehicle", massRateModelSettings, bodyMap, accelerationModelMap );

    // Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( "Vehicle" );

    Eigen::VectorXd initialBodyMasses = Eigen::VectorXd( 1 );
    initialBodyMasses( 0 ) = Vehicle.TotalMass;

    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
        std::make_shared< MassPropagatorSettings< double > >(
            bodiesWithMassToPropagate, massRateModels, initialBodyMasses, terminationSettings );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                          basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 1 ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          lvlh_to_inertial_frame_rotation_dependent_variable, "Vehicle", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(mach_number_dependent_variable, "Vehicle" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(total_acceleration_dependent_variable, "Vehicle" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(altitude_dependent_variable,"Vehicle","Earth"));

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );


    // Create propagation settings for mass and translational dynamics concurrently
    std::shared_ptr< PropagatorSettings< > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    // Define integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, fixedStepSize);
    std::cout << "Propagation properties set"<<std::endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////             PROPAGATE            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(bodyMap, integratorSettings, propagatorSettings );
    std::cout << "Propagation Complete"<<std::endl;
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputSubFolder = "Results/";

    Eigen::VectorXd finalIntegratedState = ( --integrationResult.end( ) )->second;
    // Print the position (in km) and the velocity (in km/s) at t = 0.
    std::cout << "Propagation Information" << std::endl <<
                 "The initial position vector of the vehicle is [km]:" << std::endl <<
                 systemInitialState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "The initial velocity vector of the vehicle is [km/s]:" << std::endl <<
                 systemInitialState.segment( 3, 3 ) / 1E3 << std::endl;

    // Print the position (in km) and the velocity (in km/s) at t = 100.
    std::cout << "After " << 100 <<
                 " seconds, the position vector of the vehicle is [km]:" << std::endl <<
                 finalIntegratedState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "And the velocity vector of the vehicle is [km/s]:" << std::endl <<
                 finalIntegratedState.segment( 3, 3 ) / 1E3 << std::endl;

    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "DependentVariableOutput.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "PropagationOutput.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
