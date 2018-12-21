package forDebug;

import test.AGCT;
import test.Domain;
import test.JInformationFrame;

import javax.swing.*;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;

public class Debug {

    public static void debug(Object... os) {
        String s = Arrays.deepToString(os);
        StringBuilder text = new StringBuilder();
        for (int i = 0; i < s.length(); ) {
            int j = 0;
            StringBuilder b = new StringBuilder();
            while (j < 500 && i < s.length()) {
                b.append(s.charAt(i));
                j++;
                i++;
            }
            if (i < s.length()) b.append("↓");
//            System.out.println(b.toString());
            text.append(b.toString()).append("\n");
        }
        System.out.println(text.toString());
        try{
        JInformationFrame.getInstance().appendText(text.toString());
        }catch (Exception ignore){}
    }

    private static String[] FILTER = {
            "Manifold"
    };
    private static HashSet<String> set = new HashSet<String>();

    public static void fdebug(String fileName, Object... os) {
        for (String f : FILTER) if (fileName.contains(f)) return;
        try {
            PrintWriter pw;
            if (!set.contains(fileName)) {
                set.add(fileName);
                pw = new PrintWriter(new File(fileName));
            } else pw = new PrintWriter(new FileWriter(new File(fileName), true));
            String res = Arrays.deepToString(os);
            pw.println(res);
            pw.flush();
        } catch (Exception ignored) {
        }
    }
}
