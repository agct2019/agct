package _;

import static java.util.Arrays.deepToString;

public class _ {
    private static void debug(Object... os) {
        System.err.println(deepToString(os));
    }

    public static void main(String[] args) {
        new _().run();
    }

    void run() {
        String[] ts = {
                "gb:AA266298 /DB_XREF=gi:1902420 /DB_XREF=mz63a12.r1 /CLONE=IMAGE:718078 /FEA=EST /CNT=5 /TID=Mm.2972.1 /TIER=ConsEnd /STK=0 /UG=Mm.2972 /UG_TITLE=ESTs",
                "M. musculus /GEN=beta-actin /DB_XREF=gb:M12481.1 /DEF=Mouse cytoplasmic beta-actin mRNA.",
                "", ""};
        for (String t : ts) {
            String[] ss = t.split("gb:");
            for (String s : ss) {
                if (s.contains(" ")) {
                    String u = s.split(" ")[0];
//					debug(u);
                    if (u.contains(".")) {
                        debug(u.split("\\.")[0]);
                    }
                }
            }
        }
    }
}
