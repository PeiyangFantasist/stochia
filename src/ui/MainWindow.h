#pragma once

#include "core/ProbabilityModel.h"

#include <QMainWindow>

class QLabel;
class QLineEdit;
class QListWidget;
class QPushButton;
class QComboBox;
class QSpinBox;
class QTextBrowser;
class PlotWidget;

class MainWindow final : public QMainWindow {
public:
    explicit MainWindow(QWidget* parent = nullptr);

protected:
    void closeEvent(QCloseEvent* event) override;

private:
    void buildUi();
    void buildMenus();
    QWidget* buildVariablePanel();
    QWidget* buildVisualizationPanel();
    QWidget* buildInformationPanel();
    void applyTheme();

    void addDistribution();
    void addTransformation();
    void editSelected();
    void removeSelected();
    void selectVariable(const QString& name);
    void refreshVariableList();
    void refreshPlot();
    void refreshInformation(const std::vector<double>* samples = nullptr);
    void simulate();
    void newProject();
    void openProject();
    bool saveProject();
    bool saveProjectAs();
    bool writeProject(const QString& path);
    bool readProject(const QString& path);
    bool maybeSave();
    void markModified(bool modified = true);
    QString variableDescription(const stochia::RandomVariable& variable) const;
    QString nextVariableName() const;
    void showStatus(const QString& text, bool error = false);
    void showHelp(const QString& topic = QStringLiteral("guide"));

    stochia::ProbabilityModel model_;
    QString currentPath_;
    QString selectedName_;
    bool modified_ = false;

    QListWidget* variableList_ = nullptr;
    QLabel* projectLabel_ = nullptr;
    PlotWidget* plotWidget_ = nullptr;
    QComboBox* viewCombo_ = nullptr;
    QSpinBox* sampleCount_ = nullptr;
    QLabel* plotTitle_ = nullptr;
    QLabel* typeBadge_ = nullptr;
    QLabel* expectationValue_ = nullptr;
    QLabel* varianceValue_ = nullptr;
    QLabel* sampleValue_ = nullptr;
    QTextBrowser* formulaView_ = nullptr;
    QTextBrowser* derivationView_ = nullptr;
    QLineEdit* transformName_ = nullptr;
    QLineEdit* transformExpression_ = nullptr;
    QLabel* statusLabel_ = nullptr;
};
