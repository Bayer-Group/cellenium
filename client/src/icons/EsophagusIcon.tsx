import Icon from './svg/esophagus.svg';

interface Icon {
  size: number;
}

const EsophagusIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="esophagus icon" />
    </div>
  );
};

export default EsophagusIcon;
