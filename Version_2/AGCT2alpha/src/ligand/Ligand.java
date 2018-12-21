package ligand;

public class Ligand {
    private String _name;
    private String _fileName;
    //	private String _group;
    private LigandGroup _group;
    private boolean _checked;
    private boolean _enabled;
    private boolean _discarded;

    //	public boolean isSelected(){
//		return isChecked() && isEnabled();
//	}
    public String getFileName() {
        return _fileName;
    }

    public void setFileName(String fileName) {
        _fileName = fileName;
    }

    public String getName() {
        return _name;
    }

    public String getSimpleName() {
        String[] ss = getName().split("\\.");
        if (ss.length <= 1) return getName();
        String res = "";
        for (int i = 1; i < ss.length; i++) {
            res += ss[i];
            if (i < ss.length - 1) res += ".";
        }
        return res;
    }

    public void setName(String name) {
        this._name = name;
    }

    public LigandGroup getGroup() {
        return _group;
    }

    public void setGroup(LigandGroup group) {
        this._group = group;
    }

    public boolean isChecked() {
        return _checked;
    }

    public void setChecked(boolean checked) {
        this._checked = checked;
    }

    public boolean isEnabled() {
        return !_discarded && _enabled;
    }

    public void setEnabled(boolean enabled) {
        this._enabled = enabled;
    }

    public void setDiscarded(boolean discarded) {
        this._discarded = discarded;
    }
}