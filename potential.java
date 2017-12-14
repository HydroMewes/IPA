package rub.hydrology.gamaext;

import java.util.Arrays;
import java.util.Collections;


import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import msi.gama.precompiler.GamlAnnotations.doc;
import msi.gama.precompiler.GamlAnnotations.operator;
import msi.gama.util.IList;
import msi.gama.util.matrix.GamaFloatMatrix;
import msi.gama.util.matrix.IMatrix;

import com.google.common.primitives.Doubles;


public class Potential_Implementation {

	
	@operator(value = "calcPsiM")
	@doc(value = "Return PsiM for a given soil and water content")
	public static double calcPsiM(Integer Soil, Double WaterContent, GamaFloatMatrix Soildata) {
		
		IList<Double> x = Soildata.getColumn(null, 0);
		double psiM = Double.NaN;
		
		if(WaterContent <= Collections.min(x)){
			PolynomialFunction[] splines = getFunction(Soil,Soildata).getPolynomials();
			PolynomialFunction first = splines[0];
			psiM = -1 * Math.pow(10.0,first.value(WaterContent));
		}
		
		if(WaterContent >= Collections.max(x)){
			PolynomialFunction[] splines = getFunction(Soil,Soildata).getPolynomials();
			PolynomialFunction last = splines[splines.length-1];
			psiM = -1 * Math.pow(10.0,last.value(WaterContent));
		}
		
		
		if(WaterContent > Collections.min(x) && WaterContent < Collections.max(x)){
			double pf = getFunction(Soil,Soildata).value(WaterContent);
			psiM = -1 * Math.pow(10.0,pf);
		}
		
		

		return psiM;
	}
	
	@operator(value="getPF")
	@doc("Returns pF value for a given soil and water content")
	public static double getPF(Integer Soil, Double WaterContent, GamaFloatMatrix Soildata) {
		double pf = getFunction(Soil,Soildata).value(WaterContent);
		return pf;
	}
	
	
	//Berechnet das gravimetrische Potential
		@operator(value = "getPsiZ")
		@doc("Returns gravimetric potential for a given soil and water content")
		public static double getPsiZ(double depth){
			double PsiZ = (depth * 98)/100;
			return PsiZ;
		}
		
		//Berechnet den hydraulischen Gradienten
		@operator(value = "getPsiH")
		@doc("Calculation of PsiH")
		public static double getPsiH(double psiM, double z) {
			double PsiH = Double.NaN;
			PsiH = psiM + getPsiZ(z);
			return PsiH;
		}
		
		//Berechnet kf für ungesättigte Verhältnisse
		@operator(value = "calcKF")
		@doc("Calculation of kf")
		public static double calcKF(double alpha, double ks, double n, double l, double psiM) {	
			double kf = Double.NaN;
			double m = 1 - 1 / n;
			double upper_term = Math.pow((1-Math.pow((alpha * Math.abs(psiM)),(n-1))*Math.pow((1+Math.pow((alpha*Math.abs(psiM)),n)),(-1*m))),2);
			double lower_term = Math.pow((1+ Math.pow((alpha*Math.abs(psiM)),n)),(m*l));
			
			kf = ks * upper_term/lower_term;
			
			return kf;
		}
		
		public static GamaFloatMatrix checkData(IList<Double> x_g, IList<Double> y_g) {
			
			double [] checked_data_x = {Double.NaN}; //
			double [] checked_data_y = {Double.NaN};
			
			double [] x = Doubles.toArray(x_g);
			double [] y = Doubles.toArray(y_g);
			
			int index = 9999; 
			
			int i = 0;
			int j = 0;
			
			for (double entry : y) {
			 if(Double.isNaN(entry)==true) {
					index = i;
					j = j + 1;
				}
				i = i +1 ;
			}
			
			if(j>0){
			checked_data_x = Arrays.copyOfRange(x, 0, index);
			checked_data_y = Arrays.copyOfRange(y, 0, index);
			}
			
			if(j==0) {
				checked_data_x = x;
				checked_data_y = y;
			}
			

			
			double[][] rm_data = {checked_data_x,checked_data_y};
			GamaFloatMatrix GFM_checkeddata = new GamaFloatMatrix(new Array2DRowRealMatrix(rm_data));
			return GFM_checkeddata;
		}
		
		public static PolynomialSplineFunction getFunction(Integer Soil,GamaFloatMatrix Soildata) {
			
			if(Soil > Soildata.numCols) {
				Soil = (Soildata.numCols-1);
				System.out.println("Chosen soil not available, switch to highest indexed soil in database");
			}
			
			IList<Double> x = Soildata.getColumn(null, 0);
			IList<Double> y =  Soildata.getColumn(null, Soil);
				
			GamaFloatMatrix GFM_checkedData = checkData(x, y);
			
			IMatrix<Double> T_checkData = GFM_checkedData.transpose(null);
			
			double [] x_c = Doubles.toArray(T_checkData.getColumn(null, 0));
			double [] y_c = Doubles.toArray(T_checkData.getColumn(null, 1));
			
			
			PolynomialSplineFunction pedotrsf = new LinearInterpolator().interpolate(x_c, y_c);
	
			
			return pedotrsf;
		}
 }
