package ligand;

import java.util.HashMap;

public class LigandGroup {
    private static HashMap<String, LigandGroup> map = new HashMap<String, LigandGroup>();

    public static LigandGroup get(String name) {
        if (map.containsKey(name)) return map.get(name);
        return new LigandGroup(name);
    }

    final String name;
    private boolean isSelected;

    private LigandGroup(String name) {
        this.name = name;
    }

    public boolean isSelected() {
        return isSelected;
    }

    public void setSelected(boolean isSelected) {
        this.isSelected = isSelected;
    }

    @Override
    public boolean equals(Object arg0) {
        return name.equals(((LigandGroup) arg0).name);
    }

    @Override
    public String toString() {
        return name;
    }
}
