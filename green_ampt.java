package rub.hydrology.gamaext;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import msi.gama.precompiler.GamlAnnotations.doc;
import msi.gama.precompiler.GamlAnnotations.operator;
import msi.gama.util.matrix.GamaFloatMatrix;



public class GA_Implementation {

	
	@operator(value = "getInf")
	@doc(value = "Calculates Infiltration from given event lengt and soil data")
	public static GamaFloatMatrix getInf(Integer Event_Length, double step, double eta, double ks, double Sqr_HpsiF_term) {
		
		Double event_step_no = (Event_Length/step);
		int no_steps = event_step_no.intValue();
		
		double[] t = new double[no_steps];
		
		
		double i_intern = 0.0;
		
		
		for(int i = 0; i < event_step_no.intValue(); ++i) {
			t[i] = i_intern;
			i_intern = i_intern + step;
		}
		

		
		
		double [] F_t =  new double [no_steps];
		double [] Inf =  new double [no_steps];
		double [] t_star =  new double [no_steps];

		double S_sqr = 2 * eta * ks * Sqr_HpsiF_term;

		for (int i=0; i < (no_steps); ++i) 
		{	
			t_star[i] = (2*Math.pow(ks, 2)*t[i])/S_sqr;
			double upper_term = t_star[i] + Math.log(1 + t_star[i] + (Math.sqrt(2*t_star[i]))/(1+ Math.sqrt(2*t_star[i])/6));
		    double lower_term = (1/(1 + t_star[i] + (Math.sqrt(2*t_star[i]))/(1+ Math.sqrt(2*t_star[i])/6)))-1;
		    F_t[i] = S_sqr/(2*ks) * (-1 -(upper_term/lower_term));
		}
		
		for (int i=0; i < (no_steps-1); ++i){
		    Inf[i] = (F_t[i+1]-F_t[i])/step;
		}
			
		GamaFloatMatrix GFM_int = new GamaFloatMatrix(new Array2DRowRealMatrix(Inf));
		
		return GFM_int;
	}
}
