package normalizationData;

import gui.NormalizationDataObserver;

import java.util.ArrayList;

public abstract class AbstractNormalizationData {
    private ArrayList<NormalizationDataObserver> observers = new ArrayList<NormalizationDataObserver>();

    public void addObserver(NormalizationDataObserver observer) {
        observers.add(observer);
    }

    public void deleteObserver(NormalizationDataObserver observer) {
        observers.remove(observer);
    }

    public void nortifyObservers() {
        for (NormalizationDataObserver observer : observers) {
            observer.update(this);
        }
    }
}
