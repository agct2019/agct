package ligand;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;


public class LigandList extends ArrayList<Ligand> {
    /**
     *
     */
    private static final long serialVersionUID = 3796430192162807914L;
    private LigandGroup[] _groupNames;

    public LigandGroup[] getGroups() {
        if (_groupNames != null) return _groupNames;
        final ArrayList<LigandGroup> groupnames = new ArrayList<LigandGroup>();
        for (final Ligand ligand : this) {
            final LigandGroup groupname = ligand.getGroup();
            if (groupname != null && !groupnames.contains(groupname)) {
                groupnames.add(groupname);
            }
        }
        Collections.sort(groupnames, new Comparator<LigandGroup>() {
            @Override
            public int compare(final LigandGroup arg0, final LigandGroup arg1) {
                // TODO Auto-generated method stub
                return arg0.name.compareTo(arg1.name);
            }
        });
        return _groupNames = groupnames.toArray(new LigandGroup[0]);
    }

    public int numberOfSelectedLigands() {
        int res = 0;
        for (final Ligand ligand : this)
            if (ligand.isChecked()) {
                res++;
            }
        return res;
    }
}
