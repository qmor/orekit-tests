package ru.vniiem.orekit;

import java.io.File;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataContext;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.DragSensitive;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.ThirdBodyAttraction;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.ICGEMFormatReader;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.models.earth.atmosphere.DTM2000;
import org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.events.NodeDetector;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.PVCoordinatesProvider;

public class TestProgramm {


	
	
	public static void main(String[] args) {

		System.setProperty( org.orekit.data.DataProvidersManager.OREKIT_DATA_PATH , TestProgramm.class.getResource("/orekit-data").getPath());
		File orekitData = new File(TestProgramm.class.getResource("/orekit-data").getPath());
		DataContext.getDefault().getDataProvidersManager().addProvider(new DirectoryCrawler(orekitData));

        GravityFieldFactory.clearPotentialCoefficientsReaders();
        GravityFieldFactory.addPotentialCoefficientsReader(new ICGEMFormatReader("eigen-6s.gfc", true));
        NormalizedSphericalHarmonicsProvider  hprov = GravityFieldFactory.getNormalizedProvider(36, 36);
        
		

		PVCoordinatesProvider sun = CelestialBodyFactory.getSun();
		Frame ITRFFrame = FramesFactory.getPZ9011(IERSConventions.IERS_2010 , false);
		OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, ITRFFrame);
		Frame inertialFrame = FramesFactory.getEME2000();
		TimeScale utc = TimeScalesFactory.getUTC();
		AbsoluteDate initialDate = new AbsoluteDate(2020, 01, 24, 02, 50, 42.630425, utc);

		PVCoordinates coordsECEF =  new PVCoordinates(new Vector3D(-4309.11566316175*1e3,5376.31103200996*1e3,1.41441272603228E-06*1e3) , new Vector3D(1.16764973146406*1e3,0.923459903868841*1e3,7.54353453883883*1e3));
		PVCoordinates coordsECI =  ITRFFrame.getTransformTo(inertialFrame, initialDate).transformPVCoordinates(coordsECEF);
		GeodeticPoint gp =earth.transform(coordsECEF.getPosition(), ITRFFrame, initialDate);
		Orbit initialOrbit  = new CartesianOrbit(coordsECI, inertialFrame, initialDate, org.orekit.utils.Constants.WGS84_EARTH_MU);

		// steps limits
		
		final double minStep  = 0.0001;
		final double maxStep  = 30;
		final double initStep = 30;

		// error control parameters (absolute and relative)
		final double positionError = 0.1;
		final double[][] tolerances = NumericalPropagator.tolerances(positionError, initialOrbit, initialOrbit.getType());

		// set up mathematical integrator
		AdaptiveStepsizeIntegrator integrator =  new DormandPrince853Integrator(minStep, maxStep, tolerances[0], tolerances[1]);
		integrator.setInitialStepSize(initStep);

		// set up space dynamics propagator
		NumericalPropagator propagator = new NumericalPropagator(integrator);
	
		SpacecraftState state = new SpacecraftState(initialOrbit,550);
		
		
		MarshallSolarActivityFutureEstimation msafe =
                new MarshallSolarActivityFutureEstimation("Sep2015F10\\.txt", MarshallSolarActivityFutureEstimation.StrengthLevel.STRONG);
		
		Atmosphere atmosphere = new DTM2000(msafe,  sun,earth);
		DragSensitive dragsens = new IsotropicDrag(1.0,-0.000254397782639841*10e3);
		DragForce df = new DragForce(atmosphere, dragsens);
		propagator.addForceModel(df);
		propagator.addForceModel((new ThirdBodyAttraction(CelestialBodyFactory.getSun())));
		propagator.addForceModel((new ThirdBodyAttraction(CelestialBodyFactory.getMoon())));
		propagator.addForceModel(new HolmesFeatherstoneAttractionModel(earth.getBodyFrame(), hprov));
		propagator.setInitialState(state);
		NodeDetector nodeDetector =  new NodeDetector(initialOrbit, ITRFFrame) 
		{
			int orbitNumber = 5969;
		@Override
		public Action eventOccurred(SpacecraftState s, boolean increasing) {
			if (increasing)
			{
			System.out.println(orbitNumber++);
			System.out.println("ИСК:"+s.getPVCoordinates());
			System.out.println("ГСК:"+s.getPVCoordinates(ITRFFrame));
			System.out.println(String.format("SMA %f I %f", s.getA() , Math.toDegrees(s.getI()) ));
			System.out.println(String.format("Period %f E %f", s.getKeplerianPeriod()/60.0 , Math.toDegrees(s.getE()) ));
			}

			return Action.CONTINUE;
		
		}};

		propagator.addEventDetector(nodeDetector);


		propagator.setMasterMode(30, new OrekitFixedStepHandler() {
			
			@Override
			public void handleStep(SpacecraftState currentState, boolean isLast) {
				//System.out.println(currentState.getPVCoordinates());
				
			}
		});
		state =  propagator.propagate(initialDate, new AbsoluteDate(2020, 01, 25, 23, 30, 00.000, utc));
		
	
	
	}

}