/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package garudaBinding;

import forDebug.Debug;
import jp.sbi.garuda.client.backend.BackendAlreadyInitializedException;
import jp.sbi.garuda.client.backend.GarudaClientBackend;
import jp.sbi.garuda.client.backend.ui.GarudaControlPanelFactory;
import jp.sbi.garuda.client.backend.ui.GarudaGlassPanel;
import jp.sbi.garuda.platform.commons.Gadget;
import jp.sbi.garuda.platform.commons.exception.NetworkException;
import jp.sbi.garuda.platform.commons.net.GarudaConnectionNotInitializedException;
import test.AGCT;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.io.File;
import java.util.List;

/**
 * main class of the function for garudaplatform.
 * This includes a GarudaClientBackend.ラッパークラスのようなもの。
 *
 * @author Takashi
 */
public class AGCTgaruda {

    private static final String SOFTWARE_ID = "4ad42a01-801f-47c7-8a74-b95a34a29cd8";
    private static final String SOFTWARE_NAME = "AGCT2.0";
    GarudaClientBackend backend;

    public AGCTgaruda(AGCTPropertyChangeListener parent) throws NetworkException {
        Debug.debug("BEGIN: AGCTGaruda # Constructor");
        Debug.debug("parent = ", parent);
        initGaruda(parent);
        AGCT agct = AGCT.getInstance();
        agct.setGlassPane(GarudaControlPanelFactory.getGarudaGlassPanel(backend));
        ((GarudaGlassPanel) agct.getGlassPane()).setKeyboardListeners();
    }

    /**
     * やっぱりparentを包含しなければいけないかもしれない。
     *
     * @throws NetworkException
     */
    public void initGaruda(AGCTPropertyChangeListener parent) throws NetworkException {
        Debug.debug("BEGIN: AGCTgaruda # initGaruda()", parent);
        try {
            Debug.debug("Create the instance of GarudaClientBackend");
            backend = new GarudaClientBackend(SOFTWARE_ID, SOFTWARE_NAME);
            backend.addGarudaChangeListener(parent);
            Debug.debug("Initialize the instance of GarudaClientBackend");
            backend.initialize();
            Debug.debug("Register the launch command to the instance of GarudaClientBackend");
            backend.setForceCloseOnDisconnect(false);
        } catch (GarudaConnectionNotInitializedException e) {
            System.err.println("Could not establish connection to Core.");
            System.err.println("error message: " + e);
        } catch (BackendAlreadyInitializedException e) {
            e.printStackTrace();
        }
        Debug.debug("END: AGCTgaruda # initGaruda()");
    }

    public List<Gadget> getCompatibleGadgetList() {
        return backend.getCompatibleGadgetList();
    }


    public void stop() {
        // TODO add stop command here
    }

    File fileSentToGaruda;
    public void sendFileToGaruda(File file) {
        Debug.debug("BEGIN: AGCTgaruda # sendFileToGaruda");
        try {
            fileSentToGaruda = file;
            backend.requestForLoadableGadgets("genelist", file);
        } catch (GarudaConnectionNotInitializedException e) {
            e.printStackTrace();
        }
        Debug.debug("END: AGCTgaruda # sendFileToGaruda");
    }

    File[] candidateSendingFiles = new File[0];

    public File[] getCandidateSendingFiles() {
        return candidateSendingFiles;
    }

    public void setCandidateSendingFiles(File[] files) {
        this.candidateSendingFiles = files.clone();
    }

    public File getFileSentToGaruda(final File[] files, JApplet owner) {
        if (files.length == 0) return null;
        for(File f:files)Debug.debug(f.getAbsolutePath());
        FileFilter filter = new FileFilter() {
            @Override
            public boolean accept(File file) {
                for (File f : files) if (file.getAbsolutePath().equals(f.getAbsolutePath())) return true;
                return false;
            }

            @Override
            public String getDescription() {
                return "clustering profiles";
            }
        };
        JFileChooser chooser = new JFileChooser(files[0]);
        chooser.addChoosableFileFilter(filter);
        chooser.setDialogTitle("choose a file sent to garuda.");
        chooser.showOpenDialog(owner);
        return chooser.getSelectedFile();
    }
}

