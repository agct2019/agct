package clustering;

import java.util.Arrays;

public class Cluster {
    private int captain;
    private int[] member;

    Cluster(int captain, int[] member) {
        this.captain = captain;
        this.member = member;
    }

    public int getCaptain() {
        return captain;
    }

    public int[] getMembers() {
        return member.clone();
    }

    @Override
    public String toString() {
        return String.format("captian = %d,member = %s\n", captain, Arrays.toString(member));
    }

    @Override
    public boolean equals(Object obj) {
        Cluster cluster = (Cluster) obj;
        if (captain != cluster.captain) return false;
        return Arrays.equals(member, cluster.member);
    }
}
