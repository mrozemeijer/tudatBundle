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

// Body Properties which are not in the body map
class BodyProperties {
   public:
      double referenceArea;   // Reference Area of the body
      double Mass;
};

//! Execute propagation of orbits of Apollo during entry.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::propagators;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;
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
    const double fixedStepSize = 1.0;

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
    BodyProperties FirstStage;
    BodyProperties Vehicle;
    BodyProperties UpperStage;

    // Create vehicle aerodynamic coefficients
    Vehicle.referenceArea = 50;
    FirstStage.referenceArea = 50;
    UpperStage.referenceArea = 50;

    Vehicle.Mass=5.0e3;
    FirstStage.Mass=5.0e3;
    UpperStage.Mass=5.0e3;

    bodyMap[ "Vehicle" ]->setConstantBodyMass( Vehicle.Mass );
    bodyMap[ "First Stage" ]->setConstantBodyMass( FirstStage.Mass );
    bodyMap[ "Upper Stage" ]->setConstantBodyMass( UpperStage.Mass );

    ///////////// SETTINGS FOR AERODYNAMICS VALUES ////////////
    // Define physical meaning of independent variables, in this case Mach number and angle of attack
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );

    // Define reference frame in which the loaded coefficients are defined.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;

    ///////////// DRAG FOR ENTIRE VEHICLE /////////////////////
    // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles1;
    forceCoefficientFiles1[ 0 ] = "C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CD.txt"; // Set drag coefficient file
    std::cout << "Loaded CD Value Vehicle"<< std::endl;
    forceCoefficientFiles1[ 2 ] = "C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CL.txt"; // Set lift coefficient file
    std::cout << "Loaded CL Value Vehicle"<<std::endl;

    // Load and parse files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings1 =
        readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles1, Vehicle.referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,areCoefficientsInNegativeAxisDirection );

    // Create and set aerodynamic coefficients
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings1, "Vehicle" ) );

    ///////////// DRAG FOR FIRST STAGE ///////////////////////
    // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles2;
    forceCoefficientFiles2[ 0 ] = "C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CD.txt"; // Set drag coefficient file
    std::cout << "Loaded CD Value First Stage" << std::endl;
    forceCoefficientFiles2[ 2 ] = "C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CL.txt"; // Set lift coefficient file
    std::cout << "Loaded CL Value First Stage" << std::endl;

    // Load and parse files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings2 =
        readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles2, FirstStage.referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,areCoefficientsInNegativeAxisDirection );

    // Create and set aerodynamic coefficients
    bodyMap[ "First Stage" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings2, "First Stage" ) );

    ///////////// DRAG FOR FIRST STAGE /////////////////////////
    // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles3;
    forceCoefficientFiles3[ 0 ] = "C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CD.txt"; // Set drag coefficient file
    std::cout << "Loaded CD Value Upper Stage" << std::endl;
    forceCoefficientFiles3[ 2 ] = "C:/tudatBundle/thesisTool/ThesisTool/Vehicle_CL.txt"; // Set lift coefficient file
    std::cout << "Loaded CL Value Upper Stage" << std::endl;

    // Load and parse files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings3 =
        readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles3, UpperStage.referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,areCoefficientsInNegativeAxisDirection );

    // Create and set aerodynamic coefficients
    bodyMap[ "Upper Stage" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings3, "Upper Stage" ) );

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
            std::make_shared< FromFileDataMapSettings< Eigen::Vector3d > >( cppFolder + "testThrustValues.txt" );

    // Define interpolator settings.
    std::shared_ptr< InterpolatorSettings > thrustInterpolatorSettings =
            std::make_shared< InterpolatorSettings >( linear_interpolator );

    // Create data interpolation settings
    std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolatorSettings =
            std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                thrustDataSettings, thrustInterpolatorSettings );

    // Define specific impulse
    double constantSpecificImpulse = 3000.0;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       thrustDataInterpolatorSettings, constantSpecificImpulse,
                                                       lvlh_thrust_frame, "Earth" ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for the vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.
    // Set Keplerian elements for vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 72130.0e3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.6;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 169.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 45.0 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 80.0 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 15.0 );

    // Convert vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements,
                                                                              earthGravitationalParameter );

    // Define propagation termination conditions (stop after 2 weeks).
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 4.0E6 );

    // Define settings for propagation of translational dynamics.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );

    // Crete mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = createMassRateModel( "Vehicle", std::make_shared< FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    // Create settings for propagating the mass of the vehicle
    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >( std::vector< std::string >{ "Vehicle" }, massRateModels,
                                                                  ( Eigen::Matrix< double, 1, 1 >( ) << Vehicle.Mass ).finished( ),
                                                                  terminationSettings );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                          basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          lvlh_to_inertial_frame_rotation_dependent_variable, "Vehicle", "Earth" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    // Create propagation settings for mass and translational dynamics concurrently
    std::shared_ptr< PropagatorSettings< > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    // Define integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 30.0 );
}
