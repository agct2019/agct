package timer;

import forDebug.Debug;

import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Timer {
    public static enum TimerKey {
        LoadInitialFile, Manifold, Triangulation, PCA, Clustering, ALL_CHI_2
    }

    private ArrayList<TimerKey> list = new ArrayList<TimerKey>();
    private Map<TimerKey, Long> startTime = new HashMap<TimerKey, Long>();
    private Map<TimerKey, Long> pastTime = new HashMap<TimerKey, Long>();

    public void start(TimerKey key) {
        long st = System.currentTimeMillis();
        startTime.put(key, st);
        list.add(key);
    }

    public long getPastTime(TimerKey key) {
        return System.currentTimeMillis() - startTime.get(key);
    }

    public void finish(TimerKey key) {
        Debug.debug("Timer.finish()", key.name(), getPastTime(key));
        pastTime.put(key, getPastTime(key));
    }

    public long getFinishTime(TimerKey key) {
        return pastTime.get(key);
    }

    private String tos(long time) {
        time /= 1000;
        long s = time % 60;
        time /= 60;
        long m = time % 60;
        time /= 60;
        long h = time;
        return String.format("%dh %02dm %02ds", (int) h, (int) m, (int) s);
    }

    public void showResult() {
        ArrayList<String> msg = new ArrayList<String>();
        for (TimerKey key : list) {
            if (pastTime.containsKey(key))
                msg.add(String.format("%s : %s", key.name(), tos(pastTime.get(key))));
        }
        JOptionPane.showMessageDialog(null, msg.toArray(new String[0]));
    }
}
