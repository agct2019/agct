import java.util.*;

class Feature implements Debuggable{

    public static int Max_Power_Of_Two = 8192;
    public static boolean [] trashed; //true iff a variable has been trashed
    public static Hashtable featureIndexToNonTrashedIndex; //maps non trashed index to the feature number
    public static boolean featuresSelected; //true iff features have been selected

    //Coefficients : NR & http://www.cfa.harvard.edu/~erosolow/idl/wavelet/get_coeffs.pro

    public static double C4[]=
    {0.4829629131445341,
     0.8365163037378079,
     0.2241438680420134,
     -0.1294095225512604};
    
    public static double C8[]=
    {0.230377813308,
     0.714846570552,
     0.630880767929,
     -0.027983769416,
     -0.187034811719,
     0.0308413818355,
     0.0328830116668,
     -0.0105974017850};

    public static double C10[]=
    {0.160102397974,
     0.603829269797,
     0.724308528437,
     0.138428145901,
     -0.242294887066,
     -0.032244869584,
     0.077571493840,
     -0.006241490212,
     -0.012580751999,
     0.003335725285};
    
    public static double C12[]=
    {0.111540743350, 
     0.494623890398, 
     0.751133908021,
     0.315250351709,
     -0.226264693965,
     -0.129766867567,
     0.097501605587,
     0.027522865530,
     -0.031582039318,
     0.000553842201,
     0.004777257511,
     -0.001077301085};

    public static double C14[]=
    {0.077852054085,
     0.396539319481,
     0.729132090846,
     0.469782287405,
     -0.143906003928,
     -0.224036184993,
     0.071309219266,
     0.080612609151,
     -0.038029936935,
     -0.016574541630,
     0.012550998556,
     0.000429577972,
     -0.001801640704,
     0.000353713799};

    public static double C16[]=
    {0.054415842243,
     0.312871590914,
     0.675630736297,
     0.585354683654,
     -0.015829105256,
     -0.284015542961,
     0.000472484573,
     0.128747426620,
     -0.017369301001,
     -0.044088253930,
     0.013981027917,
     0.008746094047,
     -0.004870352993,
     -0.000391740373,
     0.000675449406,
     -0.000117476784};

    public static double C18[]=
    {0.038077947363,
     0.243834674612,
     0.604823123690,
     0.657288078051,
     0.133197385824,
     -0.293273783279,
     -0.096840783223,
     0.148540749338,
     0.030725681479,
     -0.067632829061,
     0.000250947114,
     0.022361662123,
     -0.004723204757,
     -0.004281503682,
     0.001847646883,
     0.000230385763,
     -0.000251963188,
     0.000039347320};

    public static double C20[]=
    {0.026670057901, 
     0.188176800078, 
     0.527201188932,
     0.688459039454, 
     0.281172343661,
     -0.249846424327,
     -0.195946274377, 
     0.127369340336, 
     0.093057364604,
     -0.071394147166,
     -0.029457536822, 
     0.033212674059,
     0.003606553567,
     -0.010733175483, 
     0.001395351747,
     0.001992405295,
     -0.000685856695,
     -0.000116466855,
     0.000093588670,
     -0.000013264203};

    public static void flushAll(){
	trashed = null;
	featureIndexToNonTrashedIndex = null;
	featuresSelected = false;
    }

    public static void init_trashed(int d){
	trashed = new boolean[d];
	int i;
	for (i=0;i<d;i++)
	    trashed[i] = false;
    }

    public static int getNumberOfFeaturesAfterSelection(){
	if (featureIndexToNonTrashedIndex == null)
	    Matrix.perror("selected features not computed");
	return featureIndexToNonTrashedIndex.size();
    }

    public static int getNumberOfFeaturesBeforeSelection(AGCT ag){
	int i, nfeatures = 0;

	if (AGCT.Method_F == 0){
	    for (i=0;i<ag.myDomain.numberLigands;i++){
		if  ( (ag.myDomain.selectedLigand[i]) && (ag.myDomain.enabledLigand[i]) )
		    nfeatures++;
	    }
	}else if ( (AGCT.Method_F >= 1) && (AGCT.Method_F <= 10) ){
	    for (i=0;i<ag.myDomain.numberLigands;i++)
		if ( (ag.myDomain.selectedLigand[i]) && (ag.myDomain.enabledLigand[i]) ){
		    if (AGCT.Number_Of_Wavelet_Stamps == -1)
			nfeatures += waveletsStampsAutomatic(ag, ag.myDomain.domainTimes, i).length;
		    else
			nfeatures += waveletsStampsUserFixed(ag, ag.myDomain.domainTimes, i).length;
		    if ( (AGCT.Method_F == 1) && (AGCT.HAAR_WAVELETS_REMOVE_CONSTANT) )
			nfeatures--;
		}
	}else 
	    Matrix.perror("Gene :: Unauthorized Feature method");

	return nfeatures;
    }


    /****************************************************************************************
     * Haar & D4i Wavelets stuff
     *****/

    public static boolean powerOfTwo(int v){
	int i;

	i = 1;
	while (i < Max_Power_Of_Two){
	    if (v == i)
		return true;
	    i *= 2;
	}
	return false;
    }

    public static double[] waveletsStampsUserFixed(AGCT mAGCT, Vector dTimes, int nlig){
	/***********************************************************************************
	 * does the same as waveletsStampsAutomatic, but uses AGCT.Number_Of_Wavelet_Stamps
	 * to decide the time stamps
	 *****/

	if (AGCT.Number_Of_Wavelet_Stamps == -1)
	    Matrix.perror("Please choose a value <> -1 for the number of wavelet stamps");
	if (AGCT.Number_Of_Wavelet_Stamps <= 1)
	    Matrix.perror("Please choose a value > 1 for the number of wavelet stamps");

	double [] ret = new double [AGCT.Number_Of_Wavelet_Stamps];
	int nstamps = ( (Vector)  dTimes.elementAt(nlig) ).size() - 1;
	double minvalue = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(1)).doubleValue();
	double maxvalue = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(nstamps)).doubleValue();

	int i;
	for (i=0;i<ret.length;i++)
	    ret[i] = minvalue + ( ( (maxvalue - minvalue) / ( ret.length - 1 ) ) * i);

	return ret;
    }

    public static boolean checkAllWaveletStampsAutomatic(AGCT mAGCT, Vector dTimes){
	int i;
	for (i=0;i<mAGCT.myDomain.numberLigands;i++)
		if (mAGCT.myDomain.selectedLigand[i])
		    if (!checkWaveletStampAutomatic(mAGCT, dTimes, i))
			return false;
	return true;
    }

    public static boolean checkWaveletStampAutomatic(AGCT mAGCT, Vector dTimes, int nlig){
	//return false iff the program does not manage to decompose the time serie
	
	int nstamps = ( (Vector)  dTimes.elementAt(nlig) ).size() - 1;
	double delta = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(nstamps)).doubleValue();
	double vtest, vcur;
	boolean ok = true, found = false;
	int i = 0, nvalues = 1;

	do{
	    i = 0;
	    vtest = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(i+1)).doubleValue();
	    ok = true;
	    do{
		vcur = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(i+1)).doubleValue();
		while(vtest < vcur)
		    vtest += delta;
		if (vtest > vcur)
		    ok = false;
		else
		    i++;
	    }while( (ok == true) && (i < nstamps) );
	    if (ok == true)
		found = true;
	    else{
		delta /= 2.0;
		nvalues *= 2;
	    }
	}while( (found == false) && (delta > mAGCT.Very_Small) );

	return found;
    }

    public static double[] waveletsStampsAutomatic(AGCT mAGCT, Vector dTimes, int nlig){
	/************************************************************************************
	 * Returns a table in which each entry is a time step for wavelets
	 * the size of the table is a power of 2
	 * for example, 0.5, 1, 2, 4 -> [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
	 * this is done for ligand nlig
	 *****/
	
	double [] ret;
	int nstamps = ( (Vector)  dTimes.elementAt(nlig) ).size() - 1;
	double delta = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(nstamps)).doubleValue();
	double vtest, vcur;
	boolean ok = true, found = false;
	int i = 0, nvalues = 1;

	do{
	    i = 0;
	    vtest = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(i+1)).doubleValue();
	    ok = true;
	    do{
		vcur = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(i+1)).doubleValue();
		while(vtest < vcur)
		    vtest += delta;
		if (vtest > vcur)
		    ok = false;
		else
		    i++;
	    }while( (ok == true) && (i < nstamps) );
	    if (ok == true)
		found = true;
	    else{
		delta /= 2.0;
		nvalues *= 2;
	    }
	}while( (found == false) && (delta > mAGCT.Very_Small) );

	if (found == false)
	    Matrix.perror("the number of time stamps of " +  AGCTFileWriter.DATA_Ligand + " " + nlig + " does not yield a 2^k decomposition");

	ret = new double[nvalues];
	vcur = ((Double) ( (Vector)  dTimes.elementAt(nlig)).elementAt(1)).doubleValue();
	for (i=0;i<nvalues;i++){
	    ret[i] = vcur;
	    vcur += delta;
	}

	return ret;
    }

    public static void fillTimeSerie(Vector dGenes, Vector dTimes, double [] cvect, double [] begintime, double [] delta, int nstamps, int ig, int index, int j){
	//fills time serie cvect

	int k;
	double currenttime;

	currenttime = begintime[index];
	for (k=0;k<nstamps;k++){
	    cvect[k] = getCombinedValue(dTimes, j, currenttime, ((Gene) dGenes.elementAt(ig)).rawCoordinates[index], ((Gene) dGenes.elementAt(ig)).undeterminedCoordinates[index]);
	    currenttime += delta[index];
	}
    }

    public static double getCombinedValue(Vector dTimes, int nlig, double t, double [] raw, boolean [] undetermined){
	//fits the inter/intrapolation of time t to the values in raw of ligand nlig

	int i, left, right;
	double a, b, val, vl, vr;
	Vector timeStuff = (Vector) dTimes.elementAt(nlig);
	int sm = timeStuff.size() - 1;

	// seek right time

	i=0;
	right = -1;
	do{
	    val = ( (Double) timeStuff.elementAt(i+1) ).doubleValue();
	    if ( (val < t) || (undetermined[i] == true) )
		i++;
	    else
		right = i;
	}while( (i < sm) && (right == -1) );

	// seek left time

	i = sm - 1;
	left = -1;
	do{
	    val = ( (Double) timeStuff.elementAt(i+1) ).doubleValue();
	    if ( (val > t) || (undetermined[i] == true) )
		i--;
	    else
		left = i;
	}while( (i >= 0) && (left == -1) );

	if ( (left == -1) && (right == -1) )
	    Matrix.perror("Problems with time stamps in data");

	if (left == -1){
	    left = -1;
	    right = -1;
	    for (i=0;i<sm;i++){
		val = ( (Double) timeStuff.elementAt(i+1) ).doubleValue();
		if (val > t){
		    if (left == -1)
			left = i;
		    else if (right == -1)
			right = i;
		}
	    }
	}

	if (right == -1){
	    left = -1;
	    right = -1;
	    for (i = sm - 1; i >= 0; i--){
		val = ( (Double) timeStuff.elementAt(i+1) ).doubleValue();
		if (val < t){
		    if (right == -1)
			right = i;
		    else if (left == -1)
			left = i;
		}
	    }
	}

	if ( (left == -1) || (right == -1) )
	    Matrix.perror("Problems with left and right values for time stamps in data");

	if (right == left)
	    val = raw[left];
	else{
	    vr = ( (Double) timeStuff.elementAt(right+1) ).doubleValue();
	    vl = ( (Double) timeStuff.elementAt(left+1) ).doubleValue();

	    a = ( (raw[right] - raw[left]) / (vr - vl) );
	    b = raw[right] - a * vr;

	    val = ( a * t ) + b ;
	    //System.out.println("lvalue[" + vl + "] = " + raw[left] + "; rvalue[" + vr + "] = " + raw[right] + "; value[" + t + "] = " + val);
	}


	// compute value
	return val;
    }

    
    public static double [] getHaar(Gene gg, double [] v){
	int i, imax = v.length, j;
	double [] copyv = new double[imax];
	double [] valret = new double [imax];
	//!We do not remove the first coefficient, which does not traduce a variation of the expression

	do{
	    i=0;
	    j=0;
	    do{
		copyv[j] = (v[i] + v[i+1])/Math.sqrt(2.0);
		i+=2;
		j++;
	    }while(i<imax);
	    i=0;
	    do{
		copyv[j] = (v[i] - v[i+1])/Math.sqrt(2.0);
		i+=2;
		j++;
	    }while(i<imax);
	    imax /= 2;

	    /*if (gg.name.equals(new String("910"))){
	      for (i=0;i<v.length;i++)
	      System.out.print(copyv[i] + " ");
	      System.out.println("");
	      }*/

	    for (i=0;i<v.length;i++)
		v[i] = copyv[i];
	}while(imax>1);

	if (!AGCT.HAAR_WAVELETS_REMOVE_CONSTANT)
	    for (i=0;i<v.length;i++)
		valret[i] = v[i];
	else{
	    valret = new double [v.length - 1];
	    for (i=0;i<v.length-1;i++)
		valret[i] = v[i+1];
	}
	    
	return valret;
    }

    public static void filtD4i(double [] a, int n, int isign) {
	double C0=0.4829629131445341, C1=0.8365163037378077,
	    C2=0.2241438680420134, C3=-0.1294095225512603;
	double R00=0.603332511928053,R01=0.690895531839104,
	    R02=-0.398312997698228,R10=-0.796543516912183,R11=0.546392713959015,
	    R12=-0.258792248333818,R20=0.0375174604524466,R21=0.457327659851769,
	    R22=0.850088102549165,R23=0.223820356983114,R24=-0.129222743354319,
	    R30=0.0100372245644139,R31=0.122351043116799,R32=0.227428111655837,
	    R33=-0.836602921223654,R34=0.483012921773304,R43=0.443149049637559,
	    R44=0.767556669298114,R45=0.374955331645687,R46=0.190151418429955,
	    R47=-0.194233407427412,R53=0.231557595006790,R54=0.401069519430217,
	    R55=-0.717579999353722,R56=-0.363906959570891,R57=0.371718966535296,
	    R65=0.230389043796969,R66=0.434896997965703,R67=0.870508753349866,
	    R75=-0.539822500731772,R76=0.801422961990337,R77=-0.257512919478482;
	int nh,i,j;
	if (n < 8) return;
	double [] wksp = new double[n];

	nh = n >> 1;
	if (isign >= 0) {
	    wksp[0]  = R00*a[0]+R01*a[1]+R02*a[2];
	    wksp[nh] = R10*a[0]+R11*a[1]+R12*a[2];
	    wksp[1] = R20*a[0]+R21*a[1]+R22*a[2]+R23*a[3]+R24*a[4];
	    wksp[nh+1] = R30*a[0]+R31*a[1]+R32*a[2]+R33*a[3]+R34*a[4];
	    for (i=2,j=3;j<n-4;j+=2,i++) {
		wksp[i] = C0*a[j]+C1*a[j+1]+C2*a[j+2]+C3*a[j+3];
		wksp[i+nh] = C3*a[j]-C2*a[j+1]+C1*a[j+2]-C0*a[j+3];
	    }
	    wksp[nh-2] = R43*a[n-5]+R44*a[n-4]+R45*a[n-3]+R46*a[n-2]+R47*a[n-1];
	    wksp[n-2] = R53*a[n-5]+R54*a[n-4]+R55*a[n-3]+R56*a[n-2]+R57*a[n-1];
	    wksp[nh-1] = R65*a[n-3]+R66*a[n-2]+R67*a[n-1];
	    wksp[n-1] = R75*a[n-3]+R76*a[n-2]+R77*a[n-1];
	} else {
	    wksp[0] = R00*a[0]+R10*a[nh]+R20*a[1]+R30*a[nh+1];
	    wksp[1] = R01*a[0]+R11*a[nh]+R21*a[1]+R31*a[nh+1];
	    wksp[2] = R02*a[0]+R12*a[nh]+R22*a[1]+R32*a[nh+1];
	    if (n == 8) {
		wksp[3] = R23*a[1]+R33*a[5]+R43*a[2]+R53*a[6];
		wksp[4] = R24*a[1]+R34*a[5]+R44*a[2]+R54*a[6];
	    } else {
		wksp[3] = R23*a[1]+R33*a[nh+1]+C0*a[2]+C3*a[nh+2];
		wksp[4] = R24*a[1]+R34*a[nh+1]+C1*a[2]-C2*a[nh+2];
		wksp[n-5] = C2*a[nh-3]+C1*a[n-3]+R43*a[nh-2]+R53*a[n-2];
		wksp[n-4] = C3*a[nh-3]-C0*a[n-3]+R44*a[nh-2]+R54*a[n-2];
	    }
	    for (i=2,j=5;i<nh-3;i++) {
		wksp[j++] = C2*a[i]+C1*a[i+nh]+C0*a[i+1]+C3*a[i+nh+1];
		wksp[j++] = C3*a[i]-C0*a[i+nh]+C1*a[i+1]-C2*a[i+nh+1];
	    }
	    wksp[n-3] = R45*a[nh-2]+R55*a[n-2]+R65*a[nh-1]+R75*a[n-1];
	    wksp[n-2] = R46*a[nh-2]+R56*a[n-2]+R66*a[nh-1]+R76*a[n-1];
	    wksp[n-1] = R47*a[nh-2]+R57*a[n-2]+R67*a[nh-1]+R77*a[n-1];
	}
	for (i=0;i<n;i++) a[i]=wksp[i];
    }

    public static void condition(double [] a, int n, int isign) {
	double t0,t1,t2,t3;
	if (n < 4) return;
	if (isign >= 0) {
	    t0 = 0.324894048898962*a[0]+0.0371580151158803*a[1];
	    t1 = 1.00144540498130*a[1];
	    t2 = 1.08984305289504*a[n-2];
	    t3 = -0.800813234246437*a[n-2]+2.09629288435324*a[n-1];
	    a[0]=t0; a[1]=t1; a[n-2]=t2; a[n-1]=t3;
	} else {
	    t0 = 3.07792649138669*a[0]-0.114204567242137*a[1];
	    t1 = 0.998556681198888*a[1];
	    t2 = 0.917563310922261*a[n-2];
	    t3 = 0.350522032550918*a[n-2]+0.477032578540915*a[n-1];
	    a[0]=t0; a[1]=t1; a[n-2]=t2; a[n-1]=t3;
	}
    }

    public static double [] getDn(int num, double [] v){
	// num = 0 --> D4i
	// num = 1 --> D4
	// num = 2 --> D8
	// num = 3 --> D10
	// num = 4 --> D12
	// num = 5 --> D14
	// num = 6 --> D16
	// num = 7 --> D18
	// num = 8 --> D20

	int nn, n=v.length, isign = 1;
	double [] valret = new double[n];
	for (nn=0;nn<n;nn++)
	    valret[nn] = v[nn];
	
	if (n >= 4)
	    if (isign >= 0) {
		if (num==0)
		    condition(valret,n,1);
		for (nn=n;nn>=4;nn>>=1)
		    if (num==0)
			filtD4i(valret,nn,isign);
		    else if (num==1)
			filtDn(4,valret,nn,isign);
		    else if (num==2)
			filtDn(8,valret,nn,isign);
		    else if (num==3)
			filtDn(10,valret,nn,isign);
		    else if (num==4)
			filtDn(12,valret,nn,isign);
		    else if (num==5)
			filtDn(14,valret,nn,isign);
		    else if (num==6)
			filtDn(16,valret,nn,isign);
		    else if (num==7)
			filtDn(18,valret,nn,isign);
		    else if (num==8)
			filtDn(20,valret,nn,isign);
	    } else {
		for (nn=4;nn<=n;nn<<=1)
		    if (num==0)
			filtD4i(valret,nn,isign);
		    else if (num==1)
			filtDn(4,valret,nn,isign);
		    else if (num==2)
			filtDn(8,valret,nn,isign);
		    else if (num==3)
			filtDn(10,valret,nn,isign);
		    else if (num==4)
			filtDn(12,valret,nn,isign);
		    else if (num==5)
			filtDn(14,valret,nn,isign);
		    else if (num==6)
			filtDn(16,valret,nn,isign);
		    else if (num==7)
			filtDn(18,valret,nn,isign);
		    else if (num==8)
			filtDn(20,valret,nn,isign);
		if (num==0)
		    condition(valret,n,-1);
	    }

	return valret;
    }

    public static void filtDn(int num, double [] a, int n, int isign) {
	int i;
	int joff;
	int ioff = joff = -(num >> 1);
	// ioff = -2; joff = -num + 2;

	double [] cc = new double[num], cr = new double [num];
	int ncof = num;

	if (num == 4) for (i=0; i<num; i++) cc[i] = Feature.C4[i];
	else if (num == 8) for (i=0; i<num; i++) cc[i] = Feature.C8[i];
	else if (num == 10) for (i=0; i<num; i++) cc[i] = Feature.C10[i];
	else if (num == 12) for (i=0; i<num; i++) cc[i] = Feature.C12[i];
	else if (num == 14) for (i=0; i<num; i++) cc[i] = Feature.C14[i];
	else if (num == 16) for (i=0; i<num; i++) cc[i] = Feature.C16[i];
	else if (num == 18) for (i=0; i<num; i++) cc[i] = Feature.C18[i];
	else if (num == 20) for (i=0; i<num; i++) cc[i] = Feature.C20[i];
	else 
	    Matrix.perror("num not yet implemented in filtDn");
	double sig = -1.0;
	for (i=0; i<num; i++) {
	    cr[num-1-i]=sig*cc[i];
	    sig = -sig;
	}

	compFiltDn(a, n, isign, cc, cr, ioff, joff, ncof);
    }

    public static void compFiltDn(double [] a, int n, int isign, double [] cc, double [] cr, int ioff, int joff, int ncof){
	double ai,ai1;
	int i,ii,j,jf,jr,k,n1,ni,nj,nh,nmod;
	if (n < 4) return;
	double [] wksp = new double [n];
	nmod = ncof*n;
	n1 = n-1;
	nh = n >> 1;
	for (j=0;j<n;j++) wksp[j]=0.0;
	if (isign >= 0) {
	    for (ii=0,i=0;i<n;i+=2,ii++) {
		ni = i+1+nmod+ioff;
		nj = i+1+nmod+joff;
		for (k=0;k<ncof;k++) {
		    jf = n1 & (ni+k+1);
		    jr = n1 & (nj+k+1);
		    wksp[ii] += cc[k]*a[jf];
		    wksp[ii+nh] += cr[k]*a[jr];
		}
	    }
	} else {
	    for (ii=0,i=0;i<n;i+=2,ii++) {
		ai = a[ii];
		ai1 = a[ii+nh];
		ni = i+1+nmod+ioff;
		nj = i+1+nmod+joff;
		for (k=0;k<ncof;k++) {
		    jf = n1 & (ni+k+1);
		    jr = n1 & (nj+k+1);
		    wksp[jf] += cc[k]*ai;
		    wksp[jr] += cr[k]*ai1;
		}
	    }
	}
	for (j=0;j<n;j++) a[j] = wksp[j];
    }

    /***************************************************************************************
     * end of Wavelet stuff : now, slopes !
     *****/

    public static double slope(Vector dTimes, int nlig, double [] raw, boolean [] undetermined){
	  double xb=0, yb=0, cov=0, vx=0, nb=0, val=0;
	  int i;
	  Vector timeStuff = (Vector) dTimes.elementAt(nlig);

	  for (i=0;i<raw.length;i++){
	      if (undetermined[i] == false){
		  nb++;
		  xb+=( (Double) timeStuff.elementAt(i+1)).doubleValue();
		  yb+=raw[i];
	      }
	  }

	  if (nb <= 1.0)
	      Matrix.perror("Not enough values in " +  AGCTFileWriter.DATA_Ligand + " " + nlig + " to compute slopes");

	  xb/=nb;
	  yb/=nb;

	  for (i=0;i<raw.length;i++){
	      if (undetermined[i] == false){
		  cov += ( (( (Double) timeStuff.elementAt(i+1)).doubleValue() - xb) * (raw[i] - yb) );
		  vx += ( (( (Double) timeStuff.elementAt(i+1)).doubleValue() - xb)* (( (Double) timeStuff.elementAt(i+1)).doubleValue() - xb) );
	      }
	  }

	  cov/=nb;
	  vx/=nb;

	  val = cov/vx;
	  return val;
    }

    /*****************************************************************************************
     * Feature selection methods
     *****/

    public static void toFeatureIndexToNonTrashedIndex(){
	featureIndexToNonTrashedIndex = new Hashtable();
	int index = 0, i;
	for (i=0;i<trashed.length;i++){
	    if (!trashed[i]){
		featureIndexToNonTrashedIndex.put(new Integer(index), new Integer(i));
		index++;
	    }
	}
    }

    public static void noFeatureSelection(Domain d){
	int dim = 0, i; 
	for (i=0;i < ((Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumberBefore[0])).finalCoordinates.length;i++){
	    dim += ((Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumberBefore[0])).finalCoordinates[i].length;
	}
	init_trashed(dim);
	toFeatureIndexToNonTrashedIndex();
	featuresSelected = true;
	d.finishComputeGeneFeatures();
    }

    public static void featureSelection(final Domain d){
	d.myAGCT.timerTickOn();
	d.myAGCT.myInformationFrame.setText("Feature selection... " + AGCT.Feature_Selection_Method[d.myAGCT.myTabbedPane.mySelectionPane.getFeatureSelectionIndex()]);
	
	if (d.myAGCT.myTabbedPane.mySelectionPane.getFeatureSelectionIndex() == 1)
	    featureSelection_determination_coefficient(d);
	
	d.myAGCT.myInformationFrame.appendText("ok.");
	d.myAGCT.timerTickOff();
    }

    public static void featureSelection_determination_coefficient(Domain d){
	
	int dim = 0; 
	Matrix c = null;
	int i, j;
	double dv;
	int nfinal = AGCT.Max_Number_Of_Features;
	int varelim;
	
	for (i=0;i < ((Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumberBefore[0])).finalCoordinates.length;i++){
	    dim += ((Gene) d.domainGenes.elementAt(d.selectedGeneNumberToGeneNumberBefore[0])).finalCoordinates[i].length;
	}
	
	
	int ncurrent = dim;
	
	if (ncurrent < nfinal)
	    Matrix.perror("Number of final variables is > current number of variables");
	
	trashed = new boolean [dim];
	double [] averageRCurrent = new double [dim];
	double [] averageRStored = new double [dim];
	int [] index = new int [dim];
	
	d.myAGCT.dimFeatures = dim;
	d.myAGCT.generate_C(false);
	
	c = d.myAGCT.C;
	
	for (i=0;i<dim;i++){
	    trashed [i] = false;
	    for (j=0;j<dim;j++){
		dv = c.coordinates[i][j] * c.coordinates[i][j]; //DO NOT REPLACE C by THIS MATRIX at the end
		c.coordinates[i][j] = dv;
	    }
	}
	
	for (i=0;i<dim;i++){
	    averageRCurrent[i] = averageRStored[i] = 0.0;
	    index [i] = i;
	}
	
	for (i=0;i<dim;i++){
	    for (j=0;j<dim;j++)
		averageRStored[i] += c.coordinates [i][j];
	    averageRStored[i] /= (double) dim;
	    averageRCurrent[i] = averageRStored[i];
	}
	
	while(ncurrent > nfinal){
	    QuickSort.quicksort(averageRCurrent,index);
	    
	    varelim = index[index.length-1];
	    trashed[varelim] = true;
	    averageRStored[varelim] = - 2.0;
	    for (i=0;i<dim;i++){
		if (!trashed[i])
		    averageRStored[i] = ( (averageRStored[i] * (double) ncurrent) - c.coordinates[i][varelim] ) / ( (double) ncurrent - 1 );
		averageRCurrent[i] = averageRStored[i];
		index[i] = i;
	    }
	    ncurrent --;
	}
	
	//Gene.updateFinalCoordinates(this);
	
	averageRCurrent = averageRStored = null;
	index = null;
	d.myAGCT.C = null;
	
	toFeatureIndexToNonTrashedIndex();
	featuresSelected = true;
	d.finishComputeGeneFeatures();
    }
}
