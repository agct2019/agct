package test;

import clusteringProfile.JClusteringProfileFrame;
import forDebug.Debug;
import gene.Gene;
import gene.GeneList;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;


class JAGCTSelectionPane extends JPanel implements ActionListener, Debuggable {
    private static final String[] Prototype_Selection_Method = {"None", "K_Means_Cluster_Prototype_Aggregation"};

    private GeneList genes;

    public void setGenes(GeneList genes) {
        this.genes = genes;
    }

    JButton jValidateButton;
    private AGCT myAGCT;
    private Domain myDomain;
    String myName;
    ImageIcon myImageIcon;
    Box mainBox;
    private JLabel selectLabel;
    JComboBox featureMethod;
    private JComboBox featureSelectionMethod, prototypeConstructionMethod;
    JButton numberSelectedFeatures, numberWabletCoefficients, numberSelectedPrototypes;
    ArrayList<JAGCTCheckBox> checkBoxGene;
    ArrayList<JAGCTCheckBox> checkBoxLigand;

    ArrayList<JAGCTCheckBox> checkBoxGroup;

    JCheckBox keepPrototypes;

    static String actionComboFeature = "JAGCTSelectionPane_combo_select";
    static String actionComboFeatureSelection = "JAGCTSelectionPane_feature_selection";
    static String actionComboPrototypeSelection = "JAGCTSelectionPane_prototype_selection";
    static String actionValidate = "JAGCTSelectionPane_validate";
    static String String_NSF = "Max. #features : ", String_NSP = "Max. #prototypes : ", String_NWC = "#wavelet coeff. : ";

    boolean selectionAvailable;
    // false when no dataset is loaded (no selection of ligands or genes is
    // possible)

    int xMargin, yMargin, nfeat;

    JAGCTSelectionPane(AGCT a, Domain d, String name, ImageIcon ic) {
        checkBoxGene = new ArrayList<JAGCTCheckBox>();
        checkBoxLigand = new ArrayList<JAGCTCheckBox>();
        checkBoxGroup = new ArrayList<JAGCTCheckBox>();
        selectLabel = new JLabel();

        myName = name;
        mainBox = Box.createHorizontalBox();
        mainBox.setAlignmentX(LEFT_ALIGNMENT);

        selectionAvailable = false;

        myAGCT = a;
        myDomain = d;
        myImageIcon = ic;

        xMargin = 20;
        yMargin = 20;
        nfeat = 0;
    }

    public Rectangle getCaptureRectangle() {
        Rectangle bounds = getBounds();
        bounds.setLocation(getLocationOnScreen());
        return bounds;
    }

    public void toTrash() {
        removeAll();
        checkBoxGene = new ArrayList<JAGCTCheckBox>();
        checkBoxLigand = new ArrayList<JAGCTCheckBox>();
        checkBoxGroup = new ArrayList<JAGCTCheckBox>();
        selectLabel = new JLabel();
        mainBox = Box.createHorizontalBox();
        mainBox.setAlignmentX(LEFT_ALIGNMENT);
        selectionAvailable = false;
        xMargin = 20;
        yMargin = 20;
        nfeat = 0;
    }

    public int numberCurrentSelectableGenes() {
        int i, val = 0;
        for (i = 0; i < myDomain.getGenes().size(); i++)
            if (myDomain.getGenes().get(i).isSelected() && myDomain.getGenes().get(i).isEnabled())
                val++;
        return val;
    }

    public int numberCurrentSelectedGenes() {
        int i, val = 0;
        for (i = 0; i < myDomain.getGenes().size(); i++)
            if ((myDomain.getGenes().get(i).isSelected()) && myDomain.getGenes().get(i).isEnabled())
                val++;
        if ((AGCT.Max_Number_Of_Prototypes > -1) && (AGCT.Max_Number_Of_Prototypes < val))
            val = AGCT.Max_Number_Of_Prototypes;
        return val;
    }

    public int numberCurrentSelectedFeatures() {
        int val = Util_Feature.getNumberOfFeaturesBeforeSelection(myAGCT.data.getMyDomain().getLigands(), myAGCT.data.getMyDomain().getTimes());
        if ((AGCT.Max_Number_Of_Features > -1) && (AGCT.Max_Number_Of_Features < val))
            val = AGCT.Max_Number_Of_Features;
        return val;
    }

    /**
     * 単に情報を示すラベルを作って貼っているだけ。
     */
    public void computeLabel() {
        StringBuilder label = new StringBuilder();
        label.append("<HTML>");
        label.append("Summary : ");

        label.append("mean mode: " + myAGCT.getNormalizationData().getNormalizeMode().toString());
        label.append(" ; SDthreshold: " + myAGCT.getNormalizationData().getSDThreshold());

        label.append("<P>");
        label.append("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");
        label.append("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");
        label.append("&nbsp;&nbsp;");
        label.append("all #genes: " + myDomain.getGenes().size());
        label.append(" ; current #genes: " + numberCurrentSelectedGenes());
//		int s = 0;
//		for (Ligand ligand : myDomain.getLigands())
//			if (ligand.isChecked() && ligand.isEnabled())
//				s++;
        int numberOfTimes = Domain.getInstance().getTimes().length;
        label.append(" ; #time series: " + numberOfTimes);
        nfeat = numberCurrentSelectedFeatures();
        label.append("; #features: ");
        if (nfeat < 3)
            label.append("<FONT COLOR=FF0000>" + nfeat + "</FONT>");
        else
            label.append(nfeat);
        label.append(".<HTML>");
        selectLabel.setText(label.toString());
    }

    public void getReady() {
        int i;
        for (i = 0; i < myDomain.getLigands().size(); i++)
            ((JAGCTCheckBox) checkBoxLigand.get(i)).getReady();
        for (i = 0; i < myDomain.getGenes().size(); i++)
            ((JAGCTCheckBox) checkBoxGene.get(i)).getReady();
    }

    public void updateLFLigands() {
        int i;
        for (i = 0; i < myDomain.getLigands().size(); i++)
            ((JAGCTCheckBox) checkBoxLigand.get(i)).setEnabled(myDomain.getLigands().get(i).isEnabled());
    }

    public void updateLFGenes() {
        if (!AGCT.No_GeneSelection) {
            int i;
            for (i = 0; i < myDomain.getGenes().size(); i++)
                if (myDomain.getGenes().get(i).isEnabled())
                    ((JAGCTCheckBox) checkBoxGene.get(i)).setEnabled(true);
                else
                    ((JAGCTCheckBox) checkBoxGene.get(i)).setEnabled(false);
        }
    }

    public void setDomain(Domain d) {
        myDomain = d;
    }

    public void paintComponent(Graphics g) {
        if (selectionAvailable == false) {
            super.paintComponent(g);
            g.drawString("No selection available: no data loaded", xMargin, yMargin);
        }
    }

    public void setButtonLabel(JButton j, String s, int v) {
        j.setText(s + v);
    }

    public void displayComponents() {
        int id;
        Box superBox = getSuperBox();

        mainBox.add(superBox);
        mainBox.add(new Box.Filler(new Dimension(0, 0), new Dimension(1000, 0), new Dimension(1000, 0)));

        setLayout(new BorderLayout());
        add(new JScrollPane(mainBox), BorderLayout.CENTER);

        id = featureMethod.getSelectedIndex();
        AGCT.Method_F = id;
        myDomain.getGenes().updateEnabledGene();
        computeLabel();

        if (nfeat < 3) {
            featureMethod.setSelectedIndex(1);
//            setAndRecordFeature(1);
        }

        ControlProcess.put("displaySelectionPane", true);
    }

    private static ArrayList<JAGCTCheckBox> ligandSelectionBoxes = new ArrayList<JAGCTCheckBox>();

    public static ArrayList<JAGCTCheckBox> getLigandSelectionBoxes() {
        return ligandSelectionBoxes;
    }

    private Box getSuperBox() {
        Box generalBox, groupBox, superBox;
        Box geneSelectionBox, superLigandSelectionBox, featureBox, featureSelectionBox, prototypeSelectionBox, horizontalBoxNSF, horizontalBoxNSP;
        Box groupSelectionBox;
        // each element = a vertical Box containing the group elements selected
        String myBoxName;
        featureMethod = new JComboBox(AGCT.Feature_Method);
        featureMethod.setSelectedIndex(0);
        featureMethod.setActionCommand(JAGCTSelectionPane.actionComboFeature);
        featureMethod.addActionListener(this);

        numberWabletCoefficients = new JButton("");
        setButtonLabel(numberWabletCoefficients, JAGCTSelectionPane.String_NWC, AGCT.Number_Of_Wavelet_Stamps);
        numberWabletCoefficients.setActionCommand("number_wavelet_coefficients");
        numberWabletCoefficients.addActionListener(this);

        checkEnabledNWC();

        featureBox = Box.createVerticalBox();
        featureBox.setBorder(BorderFactory.createTitledBorder("1-Feature method"));
        featureBox.add(getHorizontalBoxNWC());

        featureSelectionMethod = new JComboBox(AGCT.Feature_Selection_Method);
        featureSelectionMethod.setSelectedIndex(0);
        featureSelectionMethod.setActionCommand(JAGCTSelectionPane.actionComboFeatureSelection);
        featureSelectionMethod.addActionListener(this);

        numberSelectedFeatures = new JButton("");
        setButtonLabel(numberSelectedFeatures, JAGCTSelectionPane.String_NSF, AGCT.Max_Number_Of_Features);
        numberSelectedFeatures.setActionCommand("max_number_of_features");
        numberSelectedFeatures.addActionListener(this);
        horizontalBoxNSF = Box.createHorizontalBox();
        horizontalBoxNSF.add(featureSelectionMethod);
        horizontalBoxNSF.add(numberSelectedFeatures);

        checkEnabledNSF();

        featureSelectionBox = Box.createVerticalBox();
        featureSelectionBox.setBorder(BorderFactory.createTitledBorder("2-Automatic feature selection method"));
        featureSelectionBox.add(horizontalBoxNSF);
        prototypeConstructionMethod = new JComboBox(Prototype_Selection_Method);
        prototypeConstructionMethod.setSelectedIndex(0);
        prototypeConstructionMethod.setActionCommand(JAGCTSelectionPane.actionComboPrototypeSelection);
        prototypeConstructionMethod.addActionListener(this);

        numberSelectedPrototypes = new JButton("");
        setButtonLabel(numberSelectedPrototypes, JAGCTSelectionPane.String_NSP, AGCT.Max_Number_Of_Prototypes);
        numberSelectedPrototypes.setActionCommand("max_number_of_prototypes");
        numberSelectedPrototypes.addActionListener(this);

        keepPrototypes = new JCheckBox("Keep");
        keepPrototypes.setToolTipText("Keep the prototypes computed in scenario");
        keepPrototypes.setEnabled(true);
        keepPrototypes.setSelected(true);

        horizontalBoxNSP = Box.createHorizontalBox();
        horizontalBoxNSP.add(prototypeConstructionMethod);
        horizontalBoxNSP.add(numberSelectedPrototypes);
        horizontalBoxNSP.add(keepPrototypes);

        checkEnabledNSP();

        prototypeSelectionBox = Box.createVerticalBox();
        prototypeSelectionBox.setBorder(BorderFactory.createTitledBorder("4-Automatic prototype (gene) construction method"));
        prototypeSelectionBox.add(horizontalBoxNSP);

        jValidateButton = new JButton("<HTML><i>Validate</i> my selection below</HTML>");
        jValidateButton.setActionCommand(JAGCTSelectionPane.actionValidate);
        jValidateButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        jValidateButton.addActionListener(this);

        geneSelectionBox = Box.createVerticalBox();
        geneSelectionBox.setBorder(BorderFactory.createTitledBorder("Genes"));
        geneSelectionBox.add(new JLabel("<HTML><i>By name</i></HTML>"));

        String ss = "";

        if (!AGCT.No_GeneSelection)
            for (int i = 0; i < myDomain.getGenes().size(); i++) {
                myBoxName = "gene " + myDomain.getGenes().get(i).getName();
                JAGCTCheckBox myBox = new JAGCTCheckBox(myBoxName, refGene, i);
                myBox.addItemListener(myBox);
                myBox.setSelected(myDomain.getGenes().get(i).isSelected());
                checkBoxGene.add(myBox);
                geneSelectionBox.add(myBox);
                if (i < AGCT.Max_Number_Of_Genes_Kept)
                    ss += myDomain.getGenes().get(i).toString();
            }

        JInformationFrame.getInstance().setText(ss);

        Box ligandSelectionBox = Box.createVerticalBox();
//		ligandSelectionBox.setLayout(new BoxLayout(ligandSelectionBox, BoxLayout.Y_AXIS));
//		ligandSelectionBox.setAlignmentX(TOP_ALIGNMENT);
        ligandSelectionBox.add(new JLabel("<HTML><i>By name</i></HTML>"));
//		ligandSelectionBox.setBorder(new LineBorder(Color.gray));
        for (int i = 0; i < myDomain.getLigands().size(); i++) {
            myBoxName = myDomain.getLigands().get(i).getName();
            JAGCTCheckBox myBox = new JAGCTCheckBox(myBoxName, refLigand, i);// initially true
            myBox.addItemListener(myBox);
//			myBox.setSelected(myDomain.getLigands().get(i).isChecked());
            checkBoxLigand.add(myBox);
            ligandSelectionBox.add(myBox, -1);
            ligandSelectionBoxes.add(myBox);
            if (i < myDomain.getLigands().size() - 1) {
                if (!myDomain.getLigands().get(i).getFileName().equals(myDomain.getLigands().get(i + 1).getFileName())) {
                    ligandSelectionBox.add(new JLabel(" "));
                }
            }
        }
        groupBox = Box.createHorizontalBox();
        groupBox.setBorder(BorderFactory.createTitledBorder("by group (AND)"));
//		groupBox.add(Box.createHorizontalStrut(15));
        groupSelectionBox = Box.createVerticalBox();
        groupSelectionBox.setAlignmentY(TOP_ALIGNMENT);
        groupSelectionBox.add(new JLabel("<HTML>" + "Group" + "</HTML>"));

        for (int j = 0; j < myDomain.getLigands().getGroups().length; j++) {
            myBoxName = myDomain.getLigands().getGroups()[j].toString();
            JAGCTCheckBox myBox = new JAGCTCheckBox(myBoxName, refGroup, j);
            myBox.addItemListener(myBox);
            myBox.setSelected(myDomain.getLigands().getGroups()[j].isSelected());

            checkBoxGroup.add(myBox);
            groupSelectionBox.add(myBox);
            if (j < myDomain.getLigands().getGroups().length - 1) {
                groupSelectionBox.add(new JLabel("OR"));
            }
        }
        groupBox.add(groupSelectionBox);
        groupBox.add(Box.createHorizontalStrut(15));
        geneSelectionBox.setAlignmentY(TOP_ALIGNMENT);
        ligandSelectionBox.setAlignmentY(TOP_ALIGNMENT);
        groupBox.setAlignmentY(TOP_ALIGNMENT);

        superLigandSelectionBox = new Box(BoxLayout.LINE_AXIS);
        superLigandSelectionBox.setAlignmentX(LEFT_ALIGNMENT);
        superLigandSelectionBox.setBorder(BorderFactory.createTitledBorder("Time series"));
        superLigandSelectionBox.add(ligandSelectionBox);
        superLigandSelectionBox.add(groupBox);

        generalBox = Box.createHorizontalBox();
        generalBox.add(geneSelectionBox);
        generalBox.add(superLigandSelectionBox);
        generalBox.setBorder(BorderFactory.createTitledBorder("3- User-fixed gene / time series / group selection"));
        superBox = Box.createVerticalBox();
        selectLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        selectLabel.setSize(jValidateButton.getSize());

        superBox.add(jValidateButton);
        superBox.add(selectLabel);
        superBox.add(featureBox);
        superBox.add(featureSelectionBox);
        superBox.add(prototypeSelectionBox);
        superBox.add(generalBox);
        return superBox;
    }

    private Box getHorizontalBoxNWC() {
        Box horizontalBoxNWC;
        horizontalBoxNWC = Box.createHorizontalBox();
        horizontalBoxNWC.add(featureMethod);
        horizontalBoxNWC.add(numberWabletCoefficients);
        return horizontalBoxNWC;
    }

    public void offCheckBoxes() {
        if ((checkBoxGene != null) && (checkBoxGene.size() > 0)) {
            for (int i = 0; i < checkBoxGene.size(); i++) {
                ((JAGCTCheckBox) checkBoxGene.get(i)).setEnabled(false);
            }
        }

        if ((checkBoxLigand != null) && (checkBoxLigand.size() > 0)) {
            for (int i = 0; i < checkBoxLigand.size(); i++) {
                ((JAGCTCheckBox) checkBoxLigand.get(i)).setEnabled(false);
            }
        }

        if ((checkBoxGroup != null) && (checkBoxGroup.size() > 0)) {
            for (int j = 0; j < checkBoxGroup.size(); j++) {
                checkBoxGroup.get(j).setEnabled(false);
            }
        }
    }

    public void checkEnabledNWC() {
        if (featureMethod.getSelectedIndex() == 0 || featureMethod.getSelectedIndex() == 11)
            numberWabletCoefficients.setEnabled(false);
        else
            numberWabletCoefficients.setEnabled(true);
    }

    public void checkEnabledNSF() {
        if (featureSelectionMethod.getSelectedIndex() == 0)
            numberSelectedFeatures.setEnabled(false);
        else
            numberSelectedFeatures.setEnabled(true);
    }

    public void checkEnabledNSP() {
        if (prototypeConstructionMethod.getSelectedIndex() == 0)
            numberSelectedPrototypes.setEnabled(false);
        else
            numberSelectedPrototypes.setEnabled(true);
    }

    public void finishValidate() {
        int i;

        ControlProcess.put("geneFeaturesComputed", true);
        myAGCT.computeFeatureNames();

        i = 0;
        while (myDomain.getGenes().get(i).isSelected() == false)
            i++;

        JInformationFrame.getInstance().setText(((Gene) myDomain.getGenes().get(i)).coordinatesString());

        JClusteringProfileFrame.fillTime_Stamps_Summaries(myDomain);
    }

    /**
     * Validateを行う。重要。
     */
    public void goValidate() {
        if (AGCT.MYDEBUG)
            AGCT.debug("JAGCTSelectionPane.goValidate()");
        //------ 単にボタン等を使えなくしているだけ。
        featureMethod.setEnabled(false);
        featureSelectionMethod.setEnabled(false);
        numberSelectedFeatures.setEnabled(false);
        numberWabletCoefficients.setEnabled(false);
        prototypeConstructionMethod.setEnabled(false);
        numberSelectedPrototypes.setEnabled(false);
        keepPrototypes.setEnabled(false);
        if (keepPrototypes.isSelected())
            Scenario.add("JAGCTSelectionPane_Keep_Prototypes");

        offCheckBoxes();

        //
        myDomain.updateFilterData();


        myDomain.fillSelectedGenesAllAnnotations();

        Thread t = new Thread() {
            public void run() {
                jValidateButton.setText("<HTML>Selection is <i>Validated</i></HTML>");
                jValidateButton.setEnabled(false);
            }
        };
        t.start();

        ControlProcess.put("filtersSelected", true);

        myDomain.featureAndPrototypeSelection();

    }

    // methods related to feature selection

    public void goComboFeatureSelection(ActionEvent e) {
        int id;
        id = ((JComboBox) e.getSource()).getSelectedIndex();
        setMethod_FS(id);
    }

    public void playMethod_FS(int v) {
        Debug.debug("JAGCTSelectionPane # playMethod_FS", v);
        featureSelectionMethod.setSelectedIndex(v);
        setMethod_FS(v);
    }

    public void setMethod_FS(int v) {
        AGCT.Method_FS = v;

        if (v == 0) {
            AGCT.Max_Number_Of_Features = -1;
            setButtonLabel(numberSelectedFeatures, JAGCTSelectionPane.String_NSF, AGCT.Max_Number_Of_Features);
        }
        checkEnabledNSF();
    }

    public int getFeatureSelectionIndex() {
        return featureSelectionMethod.getSelectedIndex();
    }

    public void moreMax_Number_Of_Features() {
        if (AGCT.Max_Number_Of_Features < 1) {
            AGCT.Max_Number_Of_Features = -1;
            featureSelectionMethod.setSelectedIndex(0);
        }
        setButtonLabel(numberSelectedFeatures, JAGCTSelectionPane.String_NSF, AGCT.Max_Number_Of_Features);
        computeLabel();
    }

    // methods related to prototype selection

    public void goComboPrototypeSelection(ActionEvent e) {
        int id = ((JComboBox) e.getSource()).getSelectedIndex();
        setMethod_PS(id);
    }

    public void playMethod_PS(int v) {
        prototypeConstructionMethod.setSelectedIndex(v);
        setMethod_PS(v);
    }

    public void setMethod_PS(int v) {
        AGCT.Method_PS = v;

        if (v == 0) {
            AGCT.Max_Number_Of_Prototypes = -1;
            setButtonLabel(numberSelectedPrototypes, JAGCTSelectionPane.String_NSP, AGCT.Max_Number_Of_Prototypes);
        }

        checkEnabledNSP();
    }

    public int getPrototypeSelectionIndex() {
        return prototypeConstructionMethod.getSelectedIndex();
    }

    public void moreMax_Number_Of_Prototypes() {
        if (AGCT.Max_Number_Of_Prototypes < 1) {
            AGCT.Max_Number_Of_Prototypes = -1;
            prototypeConstructionMethod.setSelectedIndex(0);
        }

        setButtonLabel(numberSelectedPrototypes, JAGCTSelectionPane.String_NSP, AGCT.Max_Number_Of_Prototypes);
        computeLabel();
    }

    // methods related to number of features

    public void moreNumber_Of_Wavelet_Stamps() {
        setButtonLabel(numberWabletCoefficients, JAGCTSelectionPane.String_NWC, AGCT.Number_Of_Wavelet_Stamps);
        computeLabel();
    }

    public void moreMethod_F() {
        if ((!Util_Feature.checkAllWaveletStampsAutomatic(myAGCT)) && (AGCT.Number_Of_Wavelet_Stamps == -1)) {
            myAGCT.defaultNWS();
            setButtonLabel(numberWabletCoefficients, JAGCTSelectionPane.String_NWC, AGCT.Number_Of_Wavelet_Stamps);
        }

        checkEnabledNWC();
    }


    public void setMethod_F(int id) {
        moreMethod_F();
        myDomain.getGenes().updateEnabledGene();
        updateLFGenes();
        computeLabel();
    }

    void comboFeatureEvent(ActionEvent e) {
        int id = ((JComboBox) e.getSource()).getSelectedIndex();
        setAndRecordFeature(id);
    }

    private void setAndRecordFeature(int id) {
        Debug.debug("JAGCTSelectionPane # setAndRecordFeature", id);
        AGCT.Method_F = id;
        if ((!Util_Feature.checkAllWaveletStampsAutomatic(myAGCT)) && (AGCT.Number_Of_Wavelet_Stamps == -1))
            Scenario.add("AGCT_Modification_Number_Of_Wavelet_Stamps", AGCT.DEFAULT_Number_Of_Wavelet_Stamps);
        setMethod_F(id);

        Scenario.add("JAGCTSelectionPane_Method_F", AGCT.Method_F + "");
    }

    public void actionPerformed(ActionEvent e) {
        String command = e.getActionCommand();
        if (command.equals(JAGCTSelectionPane.actionValidate)) {
            if (nfeat >= 3) {
                goValidate();
                Scenario.add("JAGCTSelectionPane_Validate", "");
            }
        } else if (command.equals(JAGCTSelectionPane.actionComboFeature)) {
            comboFeatureEvent(e);
        } else if (command.equals(JAGCTSelectionPane.actionComboFeatureSelection)) {
            goComboFeatureSelection(e);
            Scenario.add("JAGCTSelectionPane_Feature_Selection_Method", AGCT.Method_FS + "");
            if (AGCT.Method_FS == 0)
                Scenario.add("JAGCTSelectionPane_Max_Number_Of_Features", AGCT.Max_Number_Of_Features + "");
        } else if (command.equals(JAGCTSelectionPane.actionComboPrototypeSelection)) {
            goComboPrototypeSelection(e);
            Scenario.add("JAGCTSelectionPane_Prototype_Selection_Method", AGCT.Method_PS + "");
            if (AGCT.Method_PS == 0)
                Scenario.add("JAGCTSelectionPane_Max_Number_Of_Prototypes", AGCT.Max_Number_Of_Prototypes + "");
        } else if (command.equals("max_number_of_features")) {
            myAGCT.requestModificationMax_Number_Of_Features();
            moreMax_Number_Of_Features();
        } else if (command.equals("max_number_of_prototypes")) {
            myAGCT.requestModificationMax_Number_Of_Prototypes();
            moreMax_Number_Of_Prototypes();
        } else if (command.equals("number_wavelet_coefficients")) {
            myAGCT.requestModificationNumber_Of_Wavelet_Stamps();
            moreNumber_Of_Wavelet_Stamps();
        }
    }

    static String refGene = "JAGCTCheckBox_gene";
    static String refGroup = "JAGCTCheckBox_group";
    static String refLigand = "JAGCTCheckBox_ligand";

    class JAGCTCheckBox extends JCheckBox implements ItemListener {
        String publicName;
        String myName;
        int myI;
        private boolean ready = true;

        JAGCTCheckBox(String pn, String mn, int i) {
            super(pn, true);
            myName = mn;
            myI = i;
        }

        public void getReady() {
            ready = true;
        }

        public void itemStateChanged(ItemEvent e) {
            int iSel;

            if (ready) {
//				if (AGCT.MYDEBUG)
//					AGCT.debug(myName);
                if (myName.equals(refGene)) {
                    myDomain.getGenes().get(myI).setSelected(isSelected());
                } else if (myName.equals(refGroup)) {
                    myDomain.getLigands().getGroups()[myI].setSelected(isSelected());
                    myDomain.updateEnabledLigand();
                    myDomain.getGenes().updateEnabledGene();
                    updateLFLigands();
                    updateLFGenes();
                } else if (myName.equals(refLigand)) {
                    Scenario.add("JAGCTClickLigandCheckBox", myI);
                    myDomain.getLigands().get(myI).setChecked(isSelected());
                    myDomain.getGenes().updateEnabledGene();
                    updateLFGenes();
                } else
                    Matrix.perror("JCheckBox's name is neither gene nor group\n");
            }
            computeLabel();

            if (isSelected())
                iSel = 1;
            else
                iSel = 0;
            // TODO bugFix
            if (false)
                if (ready)
                    Scenario.add("JAGCTCheckBox_StateChanged", myName + Scenario.myToken + myI + Scenario.myToken + iSel);
        }
    }
}