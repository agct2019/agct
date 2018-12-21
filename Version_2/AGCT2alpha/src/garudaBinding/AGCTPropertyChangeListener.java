package garudaBinding;

import forDebug.Debug;
import jp.sbi.garuda.client.backend.GarudaClientBackend;
import jp.sbi.garuda.client.backend.listeners.GarudaBackendPropertyChangeEvent;
import jp.sbi.garuda.client.backend.ui.GarudaControlPanel;
import jp.sbi.garuda.client.backend.ui.GarudaGlassPanel;
import jp.sbi.garuda.platform.commons.Gadget;
import jp.sbi.garuda.platform.commons.exception.NetworkException;
import jp.sbi.garuda.platform.gadget.CoreClientAPI;
import test.AGCT;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.List;

public class AGCTPropertyChangeListener implements PropertyChangeListener {

    private static final String STATUS_LABEL_HEADER = "Status: ";
    public AGCTgaruda agctbackend;
    private AGCT owner;

    /**
     * @param owner - the owner pane
     */
    public AGCTPropertyChangeListener(AGCT owner) {
        Debug.debug("BEGIN: AGCTPropertyChangeListener # Constructor");
        try {
            agctbackend = new AGCTgaruda(this);
        } catch (NetworkException n) {
            Debug.debug("cannot start garuda", n);
            AGCT.debug("cannot start garuda");
        }
        this.owner = owner;
        Debug.debug("END: AGCTPropertyChangeListener # Constructor");
    }

    @Override
    public void propertyChange(PropertyChangeEvent pce) {
        GarudaBackendPropertyChangeEvent garudaPropertyEvt = (GarudaBackendPropertyChangeEvent) pce;

        String propertyName = garudaPropertyEvt.getPropertyName();
        if (propertyName.equals(GarudaClientBackend.LOAD_DATA_PROPERTY_CHANGE_ID)) {
            addFileContentsToInput((String) garudaPropertyEvt.getSecondProperty());
        } else if (propertyName.equals(GarudaClientBackend.LOAD_GADGET_PROPERTY_CHANGE_ID)) {
            Gadget loadableGadget = (Gadget) garudaPropertyEvt.getFirstProperty();
            Debug.debug(STATUS_LABEL_HEADER + "Loaded Gadget \"" + loadableGadget.getName() + "\"");
            String launchPath = (String) garudaPropertyEvt.getSecondProperty();
            Debug.debug("AGCT just received a message from a gadget:\n" + launchPath);
            CoreClientAPI.getInstance().sentLoadGadgetResponseToCore(true, loadableGadget);
        } else if (propertyName.equals(GarudaClientBackend.CONNECTION_NOT_INITIALIZED_ID)) {
            this.agctbackend = null;
        } else if (propertyName.equals(GarudaClientBackend.GOT_GADGETS_PROPERTY_CHANGE_ID)) {
            System.out.println(agctbackend.getCompatibleGadgetList().size());
            ((GarudaGlassPanel) owner.getGlassPane()).showPanel(agctbackend.getCompatibleGadgetList());
        } else if (propertyName.equals(GarudaControlPanel.GARUDA_CONTROL_PANEL_GADGET_SELECTED)) {
            File fileToSent = agctbackend.fileSentToGaruda;
            if (fileToSent != null) {
                ((GarudaGlassPanel) owner.getGlassPane()).hidePanel();
                //noinspection unchecked
                for (Gadget gadget : (List<Gadget>) garudaPropertyEvt.getFirstProperty()) {
                    try {
                        agctbackend.sendFileToGaruda(fileToSent);
                        agctbackend.backend.sentFileToGadget(gadget, fileToSent);
                    } catch (IllegalStateException e) {
                        System.err.println("Backend not Initialized");
                    } catch (NetworkException e) {
                        System.err.println("Connection with Garuda lost. Please make sure Garuda is running and reconnect.");
                    }
                }
            }
        }
    }

    /**
     * これからこのメソッドを改造して、ロード関数の代わりにする。
     *
     * @param filePath
     */
    private void addFileContentsToInput(String filePath) {
        owner.goLoadDomain(filePath);
    }
}
